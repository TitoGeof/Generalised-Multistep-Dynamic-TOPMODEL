function alpha=subsurface_hydraulic_gradient(DEM,FDs,beta,cs,Href)
%DEM cell IDs that aren't NaNs 
ID=find(not(isnan(DEM.Z)));
% creat a sparse DEM for speed
dem0=DEM.Z;
dem0(isnan(dem0))=0;
dem0=sparse(dem0);
%temporary storage of subsurface gradient
temp=0*ID;
%loop through cell IDs
% parfor ii=1:length(ID)
parfor ii=1:length(ID)
    temp(ii,1)=tempFun(ID(ii),cs,Href,dem0,FDs);
end
alpha=beta;
alpha(ID)=temp;
%**************************************************************************
function temp=tempFun(ID,cs,Href,dem0, FDs)
    %obtain downslope flow path for each cell ID
    [ixchannel] = flowpathextract(FDs,ID);
    %obtain elevation along the path
    H=dem0(ixchannel);
    %obtain elevation gradient (representative of hydraulic gradient)
    dH=H(1)-H;
    temp=atand(dH(1)/cs);
    if nnz(dH>1)
        dH=[dH(2:end);dH(end)];
        dH=full(dH);
        %number of DEM cells in the path
        N=length(dH(:));
        %length along the path: number of cells times cell size
        L=(1:N)'*cs;
        %refine by interpolating
        Lf=linspace(min(L(:)),max(L(:)),10*N);
        dHf=interp1(L,dH,Lf,'linear');
        %check how many cells along the path the 'd' drop in gradient is reached
        Nid=find(dHf>=Href);
        %make sure it's not empty
        if nnz(Nid)>1 
            temp=atand( dHf(Nid(1))./Lf(Nid(1)) );
        end
    end
%**************************************************************************