# The ARVE-DGVM code archive

This is the original version of the ARVE Dynamic Global Vegetation Model (ARVE-DGVM) that was developed in 2006-2012 by Jed O. Kaplan and Joe R. Melton. During this time, Kaplan and Melton working at the University of Bern, the University of Victoria, and the EPFL. The ARVE-DGVM was intended to have intermediate complexity somewhere between the LPJ-DGVM (Sitch et al., 2003) and the Community Land Model (CLM, versions 3.0, 3.5, and 4.0).

The ARVE-DGVM did work, but was computationally relatively intensive and provided little performance increase over CLM. Further development of ARVE-DGVM largely stopped when the developers moved on to new jobs and other research priorities.

One of the goals of ARVE-DGVM was to have a flexible day-night timestep for calculating all biogeophysical and biogeochemical processes, and a daily timestep for allocation and other vegetation dynamics. The day-night timestep was to facilitate comparison between the model results and carbon and water fluxes measured with eddy covariance, which are typically more reliable during daytime (when there is more turbulence and vertical transport) than at night.

This version of ARVE-DGVM is archived on GitHub for archival and research purposes. It is not supplied with input data or even a Makefile, as the model is not intended to be run. Most of the code modules were last updated between July 2010 and June 2011, which was the last period of active development on the code.

Jed Kaplan (The University of Hong Kong) and Joe Melton (Environment Canada)

Last update: April 2021
