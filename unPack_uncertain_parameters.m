%**************************************************************************
function [d,Tmax,ep,Smax,mannNhs,mannNch,Hmax,ABSTOL]=unPack_uncertain_parameters(params)
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
%average max subsurface storage [m]
Hmax    = params(7);
%--------------------------------------------------------------------------
%                 define appropreiate ode-solver tolerance values
%--------------------------------------------------------------------------
%define "AbsTol" for the ode solver
MIN         = min([Tmax;ep;mannNhs;mannNch;Hmax]);
ABSTOL      = max(min(MIN,1e-6),1e-16);
