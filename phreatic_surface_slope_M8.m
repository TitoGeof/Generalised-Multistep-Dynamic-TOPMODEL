function [alpha,e,ic,icd] = phreatic_surface_slope_M8(DEM,cs,Href)
%nan IDs
nanID     = isnan(DEM.Z);
%non-NaN cell ids to be evaluated
ID        = find(not(nanID));
%evelation data matrix
Z         = DEM.Z(:);
%(X,Y) coordiantes of each cell for calculating flow direction
sy        = DEM.size(1);
sx        = DEM.size(2);
y         = ((1:sy)')*cs;
x         = (1:sx)*cs;
%coordinates grid
[X,Y]     = meshgrid(x,y);
%set nans
X(nanID)  = NaN;
Y(nanID)  = NaN;
%squeeze
X         = X(:);
Y         = Y(:);
%vector of 'giver' cell linear indecies
ic        = [];
%vector of 'receiver' cell linear indecies
icd       = [];
%gradient in each direction
e         = [];
%mean slope angle (in all directions)
alpha     = nan(size(DEM.Z));
alpha0    = nan(length(ID),1);
%loop through cells
N         = numel(ID);
parfor ii = 1:N
    ID0             = ID(ii);
    %calclate the angle between cells and the target cell
    dY              = Y-Y(ID0);
    dX              = X-X(ID0);
    theta           = atan2d(dY,dX);
    %classify the rays
    CL              = classifyDIR2(theta);
    %distance from target cell
    D               = sqrt(dY.^2+dX.^2);
    %elevation difference relative to target cell
    dZ              = Z(ID0)-Z;
    %initialise
    ic0             = [];
    icd0            = [];
    e0              = [];
    Href0           = 2*Href;
    %while loop to ensure always finding a downslope direction to give
    %only kicks in if Href is too large for a cell to find a distance
    %within DEM bounds that has >=Href elevation drop
    while isempty(icd0)
        %new reference head
        Href0       = Href0/2;
        %see if condition is met
        cond        = dZ>=Href0;
        %loop through 8 possible directions
        for kk=1:8
            %classify into 8 directions
            temp    = find(CL{kk} & cond);
            %ensure value is found
            [d,ix]  = min(D(temp));
            ix      = temp(ix);
             if not(isempty(ix)) && not(isnan(d))
                %update giver cellls
                ic0       = [ic0;ID0];
                %update receiver cells
                neighbour = [ID0+1 ID0+sy+1 ID0+sy ID0+sy-1 ID0-1 ID0-sy-1 ID0-sy ID0-sy+1];
                icd0      = [icd0;neighbour(kk)];
                %respective gradient of giver to receiver
                e0        = [e0;dZ(ix)./(D(ix)+eps)];
                w0        = e0./sum(e0);
            end
        end
    end
    %average gradient in all driections
    alpha0(ii)    = atand(sum(w0.*e0));
    ic            = [ic;ic0];
    icd           = [icd;icd0];
    e             = [e;e0];
end
cond=ic>numel(Z) | icd>numel(Z);
ic(cond)=[];
icd(cond)=[];
e(cond)=[];
alpha(ID)   = alpha0;
%in case receivers are nans
recNAN      = isnan(Z(icd));
ic(recNAN)  = [];
icd(recNAN) = [];
e(recNAN)   = [];
%**************************************************************************
function CL = classifyDIR2(theta)
%classify to 8 directions based on angle with the target cell
CL            = cell(8,1);
CL{1}         = (theta>= 67.5  & theta<112.5 );
CL{2}         = (theta>= 22.5  & theta<67.5  );
CL{3}         = (theta>=-22.5  & theta<22.5  );
CL{4}         = (theta>=-67.5  & theta<-22.5 );
CL{5}         = (theta>=-112.5 & theta<-67.5 );
CL{6}         = (theta>=-157.5 & theta<-112.5);
CL{7}         = (theta>= 157.5 | theta<-157.5);
CL{8}         = (theta>= 112.5 & theta<157.5 );