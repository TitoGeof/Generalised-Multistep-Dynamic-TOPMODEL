function demfs = fillsinks_SAL(dem)

% % check input arguments
% error(nargchk(1, 2, nargin))

% identify nans
Inan      = isnan(dem);
% set nans to -inf
dem(Inan) = -inf;
% fill depressions using imfill with an 8-neighborhood
demfs      = imfill(dem,8,'holes');
% nans in the dem are set to nan again
demfs(Inan) = nan;

