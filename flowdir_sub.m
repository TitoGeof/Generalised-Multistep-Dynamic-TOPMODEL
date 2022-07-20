function M = flowdir_sub(DEM,ic,icd,e)

% multiple and single flow direction algorithm for Digital Elevation Models
%
% Syntax
%
%     M = flowdir(DEM)
%     M = flowdir(DEM,'propertyname',propertyvalue,...)
%     [M,W] = flowdir(...)
%
% Description
% 
%     Multiple and single flowdirection algorithm that routes 
%     through flat terrain (not sinks).
%
% Input
%
%     X,Y       coordinate matrices created by meshgrid
%     dem       digital elevation model same size as X and Y
%
% Properties
%
% propertyname     propertyvalues
%
%     'type'            'multi' (default): multiple flowdirection (dinf)
%                       'single': single flow direction (d8). Flow occurs 
%                       only along the steepest descent
% 
%     'exponent'        exponent governing the relation between flow
%                       direction and slope. Default is 1.1, which means,  
%                       there is a nonlinear relation. You may want to 
%                       increase the exponent when flow direction should 
%                       rather follow a steepest descent (single) flow 
%                       direction (e.g. 5). This option is only effective 
%                       for multiple flowdirection.
%
%     'routeflats'      choose method to route over flats/plateaus.
%                       'route' uses the function routeflats
%                       'geodesic' uses routegeodesic (requires Matlab
%                       2011b or higher)
%                       'none' does not apply any routing through flats
%
% Output
%
%     M         flowdirection (sparse matrix)
%
% See also: FLOWobj
%
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 15. March, 2009



% p = inputParser;
% p.FunctionName = 'flowdir';
% addParameter(p,'exponent',1.1,@(x) isscalar(x));
% *********************************************************************
% normal multiple FLOW DIRECTION calculation
nrc = numel(DEM.Z);
% *********************************************************************
% flow direction matrix
M = sparse(ic,icd,double(e),nrc,nrc);
% ******************************************************************
M = spfun(@(x) x.^(1.1),M);
% ******************************************************************
M = spdiags(spfun(@(x) 1./x,sum(M,2)),0,nrc,nrc) * M;