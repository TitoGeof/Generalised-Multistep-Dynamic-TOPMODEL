%**************************************************************************
function [phi,Tmax,ep,Smax,mannNhs,mannNch,d]=unPack_uncertain_parameters(params)
%--------------------------------------------------------------------------
%                          uncertain (calibration) parameters
%--------------------------------------------------------------------------
%effective/drainable porosity, phi [m]
phi       = params(1);
%maximum transmissivity (at saturation) [m^2/s]
Tmax    = params(2);
%maximum daily evaporation rate, averaged across a year [m/day]
ep      = params(3);
%m/day to m/s
ep      = ep/24/60/60;
%maximum root-zone storage [m]
Smax    = params(4);
%Manning's n coefficient for hilslope
mannNhs = params(5);
%Manning's n coefficient for channel
mannNch = params(6);
%power low exponent of conductivity decay with depth, d [m]
d       = 10; %d=10 is equivalent to exponential decay
%NOTE: d is currently fixed to approximate exponential decay because otherwise
%it will be difficult to constrain phi, Tmax and d through calibration to observed
%discharge record only. To be able to constrain the three parameters additional data
%are needed, e.g., soil storage capacity or water-table, or soil moisture time-series, or

