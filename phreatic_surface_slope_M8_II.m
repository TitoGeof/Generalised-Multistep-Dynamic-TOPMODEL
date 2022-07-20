function [alpha,e,ic,icd]=phreatic_surface_slope_M8_II(DEM,cs,Href)
% cell ids to be evaluated
ID=find(not(isnan(DEM.Z)));
%initial purturbation matrix for flow distances
P0 = zeros(DEM.size(1),DEM.size(2));
%evelation data matrix
Z=DEM.Z;
% X Y coordiantes of each cell for calculating flow direction
y=flipud((1:DEM.size(1))')*cs;
x=(1:DEM.size(2))*cs;
[X,Y]=meshgrid(x,y);
X(isnan(DEM.Z))=nan;
Y(isnan(DEM.Z))=nan;
% vector of 'giver' cell linear indecies
ic=[];
%vector of 'receiver' cell linear indecies
% icd=ID+NaN;
icd=[];
%gradient in each direction
e=[];
%mean gradient (in all directions)
alpha=0*DEM.Z;
%loop through cells
for ii=1:numel(ID)
    %calclate the angle between cells in the matrix
    dY=Y-Y(ID(ii));
    dX=X-X(ID(ii));
    theta=atan2d(dY,dX);
    %classify into 8 directions
    I=classifyDIR2(theta);
    %distance from target cell
    D=sqrt(dY.^2+dX.^2);
    %elevation difference relative to target cell
    dZ=Z(ID(ii))-Z;
    %initialise
    G=[];
    Href0=2*Href;
    %while loop to ensure always finding a downslope direction to give
    while isempty(G)
        Href0=Href0/2;
        %see if condition is met
        cond=dZ>=Href0;
        %loop through 8 possible directions
        count=0;
        for kk=1:8
            temp=find(I{kk} & cond);
            [d,ix]=min(D(temp));
            ix=temp(ix);
            %ensure value is found
            if not(isempty(ix)) && not(isnan(d))
                %updategiver cells
                ic=[ic;ID(ii)];
                %update receiver cells
                icd=[icd;ix];
                %respective gradient of giver to receiver
                e=[e;dZ(ix)./(D(ix)+eps)];
                G=[G;dZ(ix)./(D(ix)+eps)];
            end
        end
    end
    %average gradient in all driections
    alpha(ID(ii))=atand(mean(G(:)));
end

%**************************************************************************
function [Ic]=classifyDIR2(I)
Ic=cell(8,1);
Ic{3} = (I>=-22.5  & I<22.5  );
Ic{2} = (I>= 22.5  & I<67.5  );
Ic{1} = (I>= 67.5  & I<112.5 );
Ic{8} = (I>= 112.5 & I<157.5 );
Ic{6} = (I>=-157.5 & I<-112.5);
Ic{5} = (I>=-112.5 & I<-67.5 );
Ic{4} = (I>=-67.5  & I<-22.5 );
Ic{7} = (I>= 157.5 & I<-157.5);
