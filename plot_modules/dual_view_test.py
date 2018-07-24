

import yt
import numpy as np
from io_modules.manipulate_fetch_gal_fields import import_fetch_gal,prepare_unigrid

def calculate_eigenvektoren(data, sizekpc, selected_field='density', verbose = False):

  n1,n2,n3   = np.shape(data[selected_field])

  x          = np.linspace(-sizekpc/2, +sizekpc/2,n1)
  y          = np.linspace(-sizekpc/2, +sizekpc/2,n2)
  z          = np.linspace(-sizekpc/2, +sizekpc/2,n3)
  X, Y, Z    = np.meshgrid(x,y,z)

  # set coordinates and weighting field
  mass_normed = data[selected_field].flatten()
  xx          = np.vstack([X.flatten(),Y.flatten(),Z.flatten()]).transpose()
  
  # set the centeer of mass
  center_of_mass = np.zeros((3))
  center_of_mass[0] = np.sum(X.flatten() * mass_normed)/np.sum(mass_normed)
  center_of_mass[1] = np.sum(Y.flatten() * mass_normed)/np.sum(mass_normed)
  center_of_mass[2] = np.sum(Z.flatten() * mass_normed)/np.sum(mass_normed)
  #
  if verbose:
    print 'center of mass',center_of_mass
  #
  # remove center of mass from frame coordinates
  xx = xx - center_of_mass[np.newaxis,:]
  
  # init intertia tensor
  inertia     = np.zeros((3,3))
  # compute diagonal elements
  for i in xrange(3):
    jj=list()
    for j in xrange(3):
      if(i!=j):
        jj.append(j)
    for j in jj:
      inertia[i,i] = inertia[i,i] + np.sum(mass_normed[:]*(xx[:,j]**2))
  # compute off-diagonal elements
  for i,j in zip([0,1,0],[1,2,2]):
    inertia[i,j] = -np.sum(mass_normed[:]*xx[:,i]*xx[:,j])
    inertia[j,i] = inertia[i,j]

  # get eigenvectors and values
  e_values, e_vectors = np.linalg.eig(inertia)

  if verbose:
    print 'Eigenvektoren',e_vectors
    print 'eigen values ',e_values

  return e_values, e_vectors



def plot_face_edge(isnap=28, selected_field='density', sizekpc=7., cutLow=1.e-5, f_out=None):

  if f_out is None:
    f_out='test_'+selected_field+'_out'+str(isnap)+'_yt_unit_plot.png')

  from mpl_toolkits.axes_grid import AxesGrid
  import matplotlib.pyplot as plt

  fig = plt.figure()

  grid = AxesGrid(fig, (0.075,0.075,0.85,0.85),
                  nrows_ncols = (1,2),
                  axes_pad = 0.05,
                  #label_mode = "L",
                  label_mode = "1",
                  share_all = True,
                  cbar_location="right",
                  cbar_mode="single",
                  cbar_size="3%",
                  cbar_pad= 0.01)
  #grid = [fig.add_subplot(211),fig.add_subplot(212)]

  data = import_fetch_gal(isnap=isnap)

  ds, dd = prepare_unigrid(data=data, 
                           add_unit= True, 
                           regionsize_kpc= sizekpc, 
                           debug=False)

  e_value,e_vectors = calculate_eigenvektoren(data=data,sizekpc=sizekpc)

  #los_vec = e_vectors[0,:] # face on
  #los_vec = e_vectors[1,:] # edge on

  los_vec = e_vectors[2,:] # face on (again)
  #up_vec  = e_vectors[1,:]
  #up_vec  = np.cross(los_vec,up_vec)

  vec_list = [ e_vectors[2,:],e_vectors[1,:]]

  for iplot,los_vec in zip(xrange(2),vec_list):

    prj = yt.OffAxisProjectionPlot(ds = ds, center = [0,0,0], normal = los_vec , fields= selected_field
      ,width       =(4, 'kpc')
      #,north_vector=up_vec
      ,weight_field='density'
      )

    prj.set_cmap(field=selected_field, cmap='inferno')

    if selected_field == 'density':
      selected_unit = 'Msun/pc**3'
    elif selected_field == 'h2density':
      selected_unit = '1/cm**3'
    else:
      raise TypeError('unit not implemented for the field')
    prj.set_unit(selected_field, selected_unit)
    prj.set_zlim(selected_field, cutLow * dd[selected_field].max().to(selected_unit), \
                 dd[selected_field].max().to(selected_unit))

    plot        = prj.plots[selected_field]
    plot.figure = fig
    plot.axes   = grid[iplot].axes
    plot.cax    = grid.cbar_axes[iplot]

    #prj.set_minorticks('all','off')
    # Finally, this actually redraws the plot.
    prj._setup_plots()

  for cax in grid.cbar_axes:
      cax.toggle_label(True)
      #cax.axis[cax.orientation].set_label(field_select)

  prj.save(f_out, mpl_kwargs={'bbox_inches':'tight'})
  print 'dump to ',f_out



if __name__ == '__main__':
  plot_face_edge()

