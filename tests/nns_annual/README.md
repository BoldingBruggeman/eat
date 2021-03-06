
The following is taken from the description of the test case from the GOTM git repository:

Here the annual simulation of the Northern Sea at 59 deg 20' N
and 1 deg 17' E during the year 1998 as discussed by Bolding et al. (2002)
is performed.

For this simulation, time series of surface slopes were extrapolated from 
observations during autumn 1998 based on four partial tides by means of 
harmonic analysis.  All necessary meteorological data are from the UK 
Meteorological Office Model in a 6-hourly temporal resolution (meteo.dat). 
For the period of the PROVESS campaign 
(Sep 7 - Nov 7, 1998, see Howarth et al., 2002; Bolding et al., 2002)
higher resolution data (half-hourly) from research vessels are available 
and pasted in an extended meteo data file (meteonns.dat).
See GOTM test case nns_seasonal for a shorter term scenario only including 
the campaign period. For calculating the resulting surface fluxes, the bulk 
formulae from Kondo75 are used here. For the evolution of the vertical
salinity profile, which is known to stabilise stratification during summer 
months, a relaxation to results obtained with a prognostic three-dimensional 
model of the North Sea (Pohlmann, 1996). By doing so, the horizontal advection,
which is the dominant process for salinity dynamics in the Northern North Sea, 
is parameterised.

The simulated temperature looks like:

![temp](./temp.png)

The simulated salinity looks like:

![salt](./salt.png)

