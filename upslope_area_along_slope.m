% %**************************************************************************
% function areaUp=upslope_area_along_slope(DEM,FDm,ACC,BETA,cs)
% % beta=beta+1;
% % beta=beta*0;
% %DEM cell IDs that aren't NaNs
% ID=find(not(isnan(DEM.Z)));
% %temporary storage of subsurface gradient
% temp=0*ID;
% %loop through cell IDs
% parfor ii=1:length(ID)
% % for ii=1:length(ID)
%     %obtain downslope flow path for each cell ID
%     [ixchannel] = flowpathextract(FDm,ID(ii),ACC);
%     %calculate upslope area along the slope
%     temp(ii,1)=sum( cs./( cosd(BETA(ixchannel)) ) );
% end
% areaUp=0*DEM.Z;
% areaUp(ID)=temp;
%**************************************************************************
function areaUp=upslope_area_along_slope(DEM,FDm,BETA,cs)
nanIDs=isnan(DEM.Z);
%DEM cell IDs that aren't NaNs
ID = find(not(nanIDs));
%temporary storage of subsurface gradient
temp = nan(size(ID));
parfor ii = 1:length(ID)
% parfor ii=1:length(ID)
    %obtain downslope flow path for each cell ID
    OUT     = dependencemap(FDm,ID(ii));
    %calculate upslope area along the slope
    A       = cs./cosd(BETA(OUT));
   temp(ii) = sum(A(:));
end
areaUp=0*DEM.Z;
areaUp(ID)=temp;