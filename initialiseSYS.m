%**************************************************************************
function [Sx,Su,Sw,Hmax]=initialiseSYS(Nc,phi,d)
%maximum subsurface storage
Hmax  = phi*d;
%--------------------------------------------------------------------------
%                            initialise variables
%--------------------------------------------------------------------------
Sx = zeros(Nc,1)+eps;
%maximum subsurface flow per unit area
Sw = zeros(Nc,1)+ Hmax;
%unsaturated zone initial storage 
Su = Sw*0+eps;
