#!python
#cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True

# boundscheck does not check for indexing error
# wraparound  A[-1] not allowed
# nonecheck   check for None for arguments
# cdivision   does not check for division by 0

import numpy as N
cimport numpy as N
cimport cython

from libc.math cimport exp

#-------------------
# type definitions -
#-------------------

DTYPEf = N.float64
DTYPEl = N.int_

ctypedef N.float64_t DTYPEf_t
ctypedef N.int_t     DTYPEl_t

# example function
cpdef sum_along_one_axis(
                         N.ndarray[DTYPEf_t, ndim=3] cube_in
                         ):
    """
      cube_in (n0,n1,n2) : input array
    """
    cdef:
        long n0   = N.shape(cube_in)[0]
        long n1   = N.shape(cube_in)[1]
        long n2   = N.shape(cube_in)[2]
        long k, i, j

        N.ndarray[DTYPEf_t, ndim=2] out       = N.zeros((n1, n2), dtype=DTYPEf)

    for i in range(n0):
        for j in range(n1):
            out[i,j] = 0.0


    for k in range(n2):
        for j in range(n1):
            for i in range(n0):
                out[i,j] = out[i,j] + cube_in[i,j,k]

    # out = N.sum(cube_in,axis=2)

    return N.array(out,dtype=DTYPEf)


# example function
cpdef grid_particle_mass(
                    long                           Nsize,
                    DTYPEf_t                       young_dt_Myr,
                    DTYPEf_t                       old_dt_Myr,
                    N.ndarray[DTYPEf_t, ndim=2] pos_id,
                    N.ndarray[DTYPEf_t, ndim=1] mass,
                    N.ndarray[DTYPEf_t, ndim=1] epoch,
                    N.ndarray[DTYPEf_t, ndim=2] vel,
                         ):
    """
      pos_id  (n_part,3) : ids of the particles (position relative to the field)
      mass    (n_part)   : masses of the particles
      epoch    (n_part)   : epoch of the particles (in age of universe when stars formed)

      sum to cube_in the contribution to the mass due to the particles at each location
    """
    cdef:
        N.ndarray[DTYPEf_t, ndim=3] mass_cube    = N.zeros((Nsize, Nsize, Nsize), dtype=DTYPEf)
        N.ndarray[DTYPEf_t, ndim=3] epoch_cube    = N.zeros((Nsize, Nsize, Nsize), dtype=DTYPEf)
        N.ndarray[DTYPEf_t, ndim=3] velx_cube    = N.zeros((Nsize, Nsize, Nsize), dtype=DTYPEf)
        N.ndarray[DTYPEf_t, ndim=3] vely_cube    = N.zeros((Nsize, Nsize, Nsize), dtype=DTYPEf)
        N.ndarray[DTYPEf_t, ndim=3] velz_cube    = N.zeros((Nsize, Nsize, Nsize), dtype=DTYPEf)
        N.ndarray[DTYPEf_t, ndim=3] young_cube    = N.zeros((Nsize, Nsize, Nsize), dtype=DTYPEf)
        N.ndarray[DTYPEf_t, ndim=3] old_cube    = N.zeros((Nsize, Nsize, Nsize), dtype=DTYPEf)
        long n_part = N.shape(pos_id)[0]
        long k, i, j, id_part
        DTYPEf_t maxAge

    assert n_part == len(mass)
    maxAge = N.max(epoch)

    for id_part in range(n_part):
        i              = long(pos_id[id_part,0])
        j              = long(pos_id[id_part,1])
        k              = long(pos_id[id_part,2])
        mass_cube[i,j,k]  = mass_cube[i,j,k]  + mass[id_part]
        # intensive quantities in the cubes are mass weighted
        epoch_cube[i,j,k] = epoch_cube[i,j,k] + mass[id_part] * epoch[id_part]
        velx_cube[i,j,k]  = velx_cube[i,j,k]  + mass[id_part] * vel[id_part, 0]
        vely_cube[i,j,k]  = vely_cube[i,j,k]  + mass[id_part] * vel[id_part, 1]
        velz_cube[i,j,k]  = velz_cube[i,j,k]  + mass[id_part] * vel[id_part, 2]

        if epoch[id_part] > 0:
            if maxAge - epoch[id_part] <= young_dt_Myr:
                young_cube[i,j,k] = young_cube[i,j,k] + mass[id_part]

            if maxAge - epoch[id_part] <= old_dt_Myr:
                old_cube[i,j,k] = old_cube[i,j,k] + mass[id_part]

    for k in range(Nsize):
        for j in range(Nsize):
            for i in range(Nsize):
                if(mass_cube[i,j,k] > 0 ):
                    epoch_cube[i,j,k] = epoch_cube[i,j,k]/mass_cube[i,j,k]
                    velx_cube[i,j,k]  = velx_cube[i,j,k] /mass_cube[i,j,k]
                    vely_cube[i,j,k]  = vely_cube[i,j,k] /mass_cube[i,j,k]
                    velz_cube[i,j,k]  = velz_cube[i,j,k] /mass_cube[i,j,k]

    return N.array(mass_cube,dtype=DTYPEf),N.array(epoch_cube,dtype=DTYPEf),N.array(young_cube,dtype=DTYPEf),N.array(old_cube,dtype=DTYPEf),N.array(velx_cube,dtype=DTYPEf),N.array(vely_cube,dtype=DTYPEf),N.array(velz_cube,dtype=DTYPEf)


# additional example, does a density estimator with a gaussian kernel
cpdef kde_cy_gauss(
                         N.ndarray[DTYPEf_t, ndim=1] y_in,  # weight
                         N.ndarray[DTYPEf_t, ndim=1] x_in,  # position
                         N.ndarray[DTYPEf_t, ndim=1] w_in,  # dispersion [dimension of position]
                         DTYPEf_t                    x_min, # lower limit of the output
                         DTYPEf_t                    x_max, # upper limit of the output
                         long                        n_out  # dimension of the output
                         ):
    """
      gaussian kernel density estimator

      e.g. to build a spectrum use
        y_in = luminosity
        x_in = velocity
        w_in = sound speed
    """
    cdef:
        long n_in = N.shape(y_in)[0]
        long i_in, i_out,n_check
        N.ndarray[DTYPEf_t, ndim=1] x_out       = N.zeros(n_out, dtype=DTYPEf)
        N.ndarray[DTYPEf_t, ndim=1] y_out       = N.zeros(n_out, dtype=DTYPEf)
        DTYPEf_t                    norm

    #-------------
    # SETUP part -
    #-------------

    sq_2_pi = N.sqrt(2.e0*2.e0*N.arcsin(1.e0))
    x_out = N.linspace(x_min,x_max,n_out)
    for i_out in range(n_out):
        y_out[i_out] = 0.0

    n_check = n_in
    for i_in in xrange(n_in):
      if(w_in[i_in]>0):
        n_check = n_check -1

    if(n_check != 0 ):
      print 'ERROR'

    #-------------------
    # CALCULATION part -
    #-------------------

    for i_in in range(n_in):

      #if(x_in[i_in] >= (x_min-w_in[i_in]) and x_in[i_in] <= (x_max+w_in[i_in])):

        norm = w_in[i_in]*sq_2_pi
        norm = y_in[i_in]/norm

        for i_out in range(n_out):
            tmp          = (x_out[i_out]-x_in[i_in])/w_in[i_in]
            #tmp          = -0.5*(tmp**2)
            tmp          = -0.5*tmp*tmp
            tmp          = exp(tmp)
            y_out[i_out] = y_out[i_out] +norm*tmp

    return N.array(x_out,dtype=DTYPEf),N.array(y_out,dtype=DTYPEf)

