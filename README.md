# mc_eor

See github commit and code doc for other progres.. lazy to update this README file.

Sept 10th 2018
- for paper... 
  1. finish analysis for Sec. 4, and update paper correspondingly
  2. fix cosmetics of plots (e.g., Fig 3; ranges shown in Fig 7,8,10)
    + add annotation of n_cut directly on Fig 3  
    + bigger axes label fonts in Fig 1.
    + add annotation of disk/satelites/sub-MCs directly on Fig 5-11
  3. re-arrange sections/text as warranted, based on comments from co-authors
  4. discussion of turbulence relating to observables (e.g., shift of SLED to higher J for higher \sigma and Mach number, but note the shift is not unique to just \sigma, will need dispersion map too, but then the latter is affected by LOS and inclination effects..)

July 16th 2018
- two really problematic bugs in yt clump find module
  1. sometimes annotate_clumps() will plot all clumps, but sometimes it won't even if you explicitly ask it to.. Work around: loop through all the clumps and plot one-by-one
  2. yt sometimes will reset the region to look for "children", such that the region at higher density it looks for can be greater than the parent region which is bogus. Work around: brute force loop through each density cut instead of using the yt find_clump(cmin, cmax, step) functionality.
  3. minor bug: plot_args in annotate_clumps() sometime takes 'color', sometime takes 'colors'. No errors is caught, it will just completely ignore the keyword argument..

July 11th 2018
- tracing cloud evolution. We have finely sampled snapshots of 0.1 Myr or 0.5 Myr, but we don't have tracer particles..
  + to track evolution of clouds, we could check by (w/ small delta_t between snapshots) 
    1. visualizing the clouds identified in each snapshots by plotting
    2. look at the cloud positions (e.g., fcoords or icoords) and velocity distributions, then calc. a pseudo upper limit on how far we expect the cloud to be displaced over that delta_t. Then, check neighboring cells to try to connect the clouds across snapshots.
- worry that we may cut off/break up cloud into multiple components because of the selected n_threshold..
  + we should probably look into how much variations do we see in the MC properties relations w/ different n_threshold (or n_cut).
    * probably interesting to plot as one figure in the paper.

July 8th 2018
- note about LoG and DoG stuff: http://www.cse.psu.edu/~rtc12/CSE486/lecture11_6pp.pdf

July 6th 2018
- even after identifying clumps/blobs, how can we tell if they are self-gravitating structure?
    + get the fields to calc KE, PE after to check they are self-gravitating.
    + check also alpha parameter (PE vs. Pressure).

July 5th 2018
- see progress in code documentation and github commit.

July 4th 2018
- see progress in code documentation and github commit.

July 3rd 2018
- continued from plotting contours over H2 density map..

July 2nd 2018
- need to think about how to flatten all the AMR hydro (levels) onto one, and then in a form that we can parse through clustering algorithms.
- take out the idea of using disperse after talking to M. Mac Low
- Maybe attempt with clumpfind (3D version) but be aware of the caveats
- read paper by Juan C. Ibáñez-Mejía about multiple clodus in potential wells (also careful about the filaments).
- note overlapping regions.

- science - investigate the dynamical properties of molecular cloud complexes at EoR in simulations.
  - look at molecular cloud complexes in models and their evolution (across diff. snapshots)
  - physical properties of interests: 
      + e.g., gas surface density, SK relation of the identified cloud complexes
      + velocity dispersion of clouds, size, SD
      + (thermal and non-thermal) pressure
      + etc
  - after ID-ing MCs, compare against obs. in nearby Dor 30? W49? 
  - note, we may care about using H2 as criterion at early times since it's needed to cool down to dense enough regions for SF, but soon as the first gen. of stars form, there are other species that can do the cooling, so H2 should not be used as the only defining criterion for SF. 
    - to confirm, why don't we do the ID using 1) H2 and 2) Pressure,density from the AMR output and compare the two?! 

~ 20 Mpc box
- used halo finder (FoF) to find main galaxy, center on main galaxy and extracted ~3 kpc around as the subregion..
    + half (baryonic) mass radius ~ 0.5 kpc.
    + r_200 DM mass radius ~ 15 kpc.
