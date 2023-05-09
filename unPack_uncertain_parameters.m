%**************************************************************************
function [d,Tmax,ep,Smax,mannNhs,mannNch,phi]=unPack_uncertain_parameters(params)
%--------------------------------------------------------------------------
%                          uncertain (calibration) parameters
%--------------------------------------------------------------------------
%powerlaw decay parameter [m]
d       = params(1);
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
%effective subsurface storage/porosity, phi [m]
phi    = params(7);

