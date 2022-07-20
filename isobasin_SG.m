function [B,NB]=isobasin_SG(DEM,FDs,As,cs,Athresh)
%define a channel network
ch       = 0*DEM.Z;
ch(As*cs^2>=Athresh)=1;
%find non-channel cell IDs
NONchID     = find(ch==0);
%cells receiving flow (only the unique values of 'receiver' cells)
cellID   = unique(FDs.ixc);
%exclude the ones that are out of the channel netwtork
cellID(ismember(cellID,NONchID))=[];
%upslope contributing area (single flow accumulation * cellsize^2)
As       = As(cellID)*cs^2;
%sort cells according to closeness to the threshold area
[~,inds] = sort(abs(As-Athresh));
%reorder cellIDs
cellID   = cellID(inds);
%master Basin ID
B        = zeros(size(DEM.Z));
%basin counter
bID      = 0;
%total number of cells
N        = length(cellID);
%loop through cells
for ii   = 1:N
    %only if the cell is not already in the mapped basin (master basin)
    if B(cellID(ii)) == 0 
        %get the dependence map
        B0 = dependencemap(FDs,cellID(ii));
        %exclude those already picked in the master basin file
        B0 = B0.*(1-sign(B));
        %calculate the area of the isobasin
        A = sum(B0(:)*cs^2);
        %if it is larger than the threshold*correction
        if A >= Athresh*0.90
            %add one to the basin IDs
            bID     = bID+1;
%             %save isobasin area
%             area(bID,1) = A;
            %update the master basin file
            B(B0>0) = bID;
        end
    end
end
%makes sure no cell has been left out
B           = DEM.Z*0+B;
% if yes, assign a new basin ID
if nnz(B==0)>0
    B       = B+1;
end
%set nans back to zero
B(isnan(B)) = 0;
%shuffle values and reset to range from 1:numel(R)
B           = shufflelabel(B,true);
%number of isobasins
NB          = max(unique(B(:)));

% % OLD CODE BASED ON STREAM ORDER
% %obtain stream order
% S               = streamorder(FDs,ch);
% %retain unique IDs
% B0              = drainagebasins(FDs,S,1) ;
% B               = B0+1;
% B(isnan(DEM.Z)) = NaN;
% %adjust basin ID to go from 1 to Nr (Nr would be number of channel reaches too)
% NB          = unique(B);
% NB(isnan(NB) | NB==0)=[];
% for ii=1:length(NB)
%     cond=B == NB(ii);
%     B(cond) = ii;
%     NB(ii) = ii;
% end
% %convert class to int32
% Y  = typecast(B(:),'int32');
% NB = typecast(NB,'int32');
% %reshape into matrix
% Y  = reshape(Y,(size(B)));
% %convert into double
% B  = double(Y);
% NB = double(NB);
% NB = max(NB);