# mc_eor

See github commit and code doc for other progres.. lazy to update this README file.

Sept 12th 2018
- more for paper, Mordecai's quick comments to be addressed
Reading the abstract, I see two claims that I think need better support:

1. "The more massive and bigger MCs in Althæa compared to the Milky Way likely result from the higher gas mass fraction, surface density, and velocity dispersion, which set the scale for fragmentation.”  It seems to me that the question of whether this actually occurs because of the effects of numerical resolution and the feedback approximation needs to be better investigated. For the resolution question, even though we don’t have a resolution study available, we can still ask the question of what the smallest objects that you are measuring are.  If they are already larger than the largest Milky Way MCs, then it is pretty clear that the simulation cannot address the question of whether there is a difference between Milky Way and simulated MC masses and radii.  The other issue here is the nature of the feedback, which likely cannot resolve the detailed interaction of the feedback with the high density regions, which could also result in fragmentation of clouds.

2. "The velocity dispersion remains 100 km s−1 even when we increase the ncut and even for the molecular substructures, likely resulting from the strong supernova and stellar feedback Althæa experienced over the multiple episodes of bursty star formation.”  Although you included external pressure in the analytical virial analysis of section 3.1, you did not include rotation, which must play some role in determining the apparent velocity dispersion. Demonstrating that it is or isn’t important seems necessary for this paper to be complete.  The idea that the supernovae could be responsible for the observed velocity dispersion seems to conflict somewhat with the fast decay of supersonic turbulence.  Perhaps we could do an analysis like that of Tamburro et al. 2009 to determine the expected equilibrium velocity dispersion for a given SN rate, and compare that to the values you are finding.

Another issue that comes up in the virial analysis is the suggestion that confinement by external pressure could be more effective than rotation in confining the MCs with very high virial number.  Since you have direct access to the turbulent and thermal pressure fields, you ought to be able to directly demonstrate or disprove the existence of such enormous confining pressures, which would place this idea on substantially stronger ground.


Sept 10th 2018
- for paper... 
  1. finish analysis for Sec. 4, and update paper correspondingly (e.g., sec 5.2 - origin of MCs and sec 5.3 on origin of CMF stuff)
  2. fix cosmetics of plots (e.g., Fig 3; ranges shown in Fig 7,8,10)
    + add annotation of n_cut directly on Fig 3  
    + bigger axes label fonts in Fig 1.
    + add annotation of disk/satelites/sub-MCs directly on Fig 5-11
    + bigger text for Fig 13
  3. re-arrange sections/text as warranted, based on comments from co-authors
  4. discussion of turbulence relating to observables (e.g., shift of SLED to higher J for higher \sigma and Mach number, but note the shift is not unique to just \sigma, will need dispersion map too, but then the latter is affected by LOS and inclination effects..)
  5. expand sec 5.3 a bit more, now it's not totally conherent.. I still need to flesh out some ideas and texts
  6. currently fig. 8 not used anywhere.. --> need to somehow, somewhere talk about temporal evolution in the main text (maybe in Discussions section?)

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
