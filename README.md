# SPB2_CT_Simulations

OpticlessPointSourceInFP:
Simulates a point source in the focal plane of the SPB2 camera. This is a very basic simulation. It does not take the optics into account, which means the point spread function is perfect. For each simulated event, a random number of photons is simulated that is randomly split between the two bi-focal spots in the camera. 
CARE can then simulate a finite point spread function by smearing the photon
positions in the camera over a 2D Gaussian.
All photons arrive at the same time
All photons have a wavelength of 300nm
The position of the event in the camera is picked at random and is in an area slightly larger than the camera itself.
The simulated events are stored in a .root file that is compatible with CARE.

CAREConfig:
CARE configuration files for SPB2

SimAnalysisMacros:
macros to run CARE and extract some useful information out of the simulated
results




