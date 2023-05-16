%**************************************************************************
function [d,Tmax,phi,ep,Smax,mannNhs,mannNch]=unPack_uncertain_parameters(params)
%--------------------------------------------------------------------------
%                          uncertain (calibration) parameters
%--------------------------------------------------------------------------
%power low exponent of conductivity decay with depth, d [m]
d        = params(1);
%maximum transmissivity (at saturation) [m^2/s]
Tmax     = params(2);
%effective/drainable porosity, phi [m]
phi      = params(3);
%maximum daily evaporation rate, averaged across a year [m/day]
ep       = params(4);
%m/day to m/s
ep       = ep/24/60/60;
%maximum root-zone storage [m]
Smax     = params(5);
%Manning's n coefficient for hilslope
mannNhs  = params(6);
%Manning's n coefficient for channel
mannNch  = params(7);



