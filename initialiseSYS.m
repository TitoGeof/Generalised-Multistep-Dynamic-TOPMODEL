%**************************************************************************
function [Sx,Su,Sw]=initialiseSYS(Nc)
%--------------------------------------------------------------------------
%                            initialise variables
%--------------------------------------------------------------------------
Sx = zeros(Nc,1)+eps;
%maximum subsurface flow per unit area
Sw = zeros(Nc,1)+eps;
%unsaturated zone initial storage 
Su = Sw*0+eps;
