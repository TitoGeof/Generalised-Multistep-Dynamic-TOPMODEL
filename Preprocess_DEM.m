function [WxmD,WxmU,WbmD,WbmU,D,SINa,SINb,COSa,COSb,areaf,AREA,Nc,Nr,HSUs...
           ,cs,DEMx,cW]=Preprocess_DEM(catchname,CHthresh,Href,Athresh,cs,ISOBASINS,DIFFUSION,TI_bins,outletW)
warning('off','all');
%--------------------------------------------------------------------------
%                                 load DEM 
%--------------------------------------------------------------------------
%load surface DEM [object]              
DEMx                 = GRIDobj([pwd '\DATA\' catchname 'dem.tif']);
%just to be sure, sometimes it comes out as 'single' and that's an issue
DEMx.Z               = double(DEMx.Z);
%fill the sinks 
DEMx                 = fillsinks(DEMx);
%turn the DEM inside out to calculate upslope adjacency matrix
DEMxi                = Invert_DEM(DEMx);
%--------------------------------------------------------------------------
%                             flow objects 
%--------------------------------------------------------------------------
%single flow direction for surface[object]
FDxs                 = FLOWobj(DEMx,'type','single');
% %multiple flow direction fo surface [object]
FDxm                 = FLOWobj(DEMx,'type','multi');
%multiple flow direction flow matrix for surface in downslope direction
MxmD                 = flowdir(DEMx,'type','multi','routeflats','route');
%single flow direction flow matrix for surface in downslope direction
MxsD                 = flowdir(DEMx,'type','multi','routeflats','route','exponent',10);
%multiple flow direction flow matrix for surface in upslope direction
MxmU                 = flowdir_inverse(DEMxi,'type','multi','routeflats','route');
%single flow direction flow matrix for surface in upslope direction
MxsU                 = flowdir_inverse(DEMxi,'type','multi','routeflats','route','exponent',10);
%single flow accumulation for surface
Axs                  = flowacc(FDxs);
%multiple flow accumulation for surface
Axm                  = flowacc(FDxm);
%surface gradient [object]
GRADx                = gradient(FDxm,DEMx);
%convert to degrees 
betaD                = atand(GRADx.Z);
%average downslope cell angles [degrees]
betaD                = betaD + 0.01;
betaD(isnan(DEMx.Z)) = NaN;
%--------------------------------------------------------------------------
%                    subsurface hydraulic gradient 
%--------------------------------------------------------------------------
if strcmp(DIFFUSION,'on')
  %NOTE: this can take a long time to run for large DEMs, so only uncomment for
  %deep systems, or flat systems, where phreatic surface slope is expected
  %to differ significantly from surface gradient (see model's journal article)
  
  %once the files are save, you can just load them if you need to run the 
  % Preprocess_DEM code for different CHthresh, or TI classification,...
  
  %reconstruct water table gradient
  [alphaD,eD,icD,icdD] = phreatic_surface_slope_M8(DEMx,cs,Href);
  save([pwd '\DATA\' catchname '_alpha_e_ic_icd_downslope' '_' num2str(Href) 'm.mat'],'alphaD','eD','icD','icdD')
  load([pwd '\DATA\' catchname '_alpha_e_ic_icd_downslope' '_' num2str(Href) 'm.mat'],'alphaD','eD','icD','icdD')
  %avoid zeros, add numbers of order 10e-2
  alphaD = alphaD + 0.01;
  %multiple flow direction flow matrix for subsurface in downslope direction
  MbmD                 = flowdir_sub(DEMx,icD,icdD,eD);
  %subaurfce gtradient in upslope direction
  [alphaU,eU,icU,icdU] = phreatic_surface_slope_M8(DEMxi,cs,Href);
  save([pwd '\DATA\' catchname '_alpha_e_ic_icd_upslope' '_' num2str(Href) 'm.mat'],'alphaU','eU','icU','icdU')
  load([pwd '\DATA\' catchname '_alpha_e_ic_icd_upslope' '_' num2str(Href) 'm.mat'],'alphaU','eU','icU','icdU')
  alphaU = alphaU + 0.01;
  %multiple flow direction flow matrix for subsurface in upslope direction
  MbmU                 = flowdir_sub_inverse(DEMxi,icD,icdD,eD);
else
  alphaD=betaD;
  alphaU=betaD;
  MbmD=MxmD;
  MbmU=MxmU;
end
%--------------------------------------------------------------------------
%                    iso-basins and channel segments
%--------------------------------------------------------------------------
if strcmp(ISOBASINS,'on')
  %once you've ran and saved is-basins, you can just load them to save time
  [B,NB]               = isobasin_SG(DEMx,FDxm,Axs,cs,Athresh);
  save([pwd '\DATA\' catchname '_isobasins_map' num2str(Athresh) '.mat'],'B','NB')
  load([pwd '\DATA\' catchname '_isobasins_map' num2str(Athresh) '.mat'],'B','NB')
else
  B=sign(DEMx.Z);
  B(isnan(B))=0;
  NB=max(B(:));
end
%channel network (ch) and channel reach IDs (R)
[R,ch]               = channel_reaches(DEMx,FDxs,CHthresh);
%--------------------------------------------------------------------------
%                 define Hydrologically Similar Units (HSUs)
%--------------------------------------------------------------------------
%upslope contributing area [m^2]
areaUp               = Axs.*cs;
%Tpographic Wetness Index (TI)
TpInd                = log(areaUp./sind(alphaD)+eps);
TpInd(isnan(DEMx.Z)) = NaN;
MAX                  = max(TpInd(:));
TpInd(betaD<=0.01)   = MAX;
%--------------------------------------------------------------------------
%                 based on Glossop channel bin edges
%--------------------------------------------------------------------------
%obtain class IDs and related slope informaiton
[HSUs,SINa,SINb,COSa,COSb,areaf,AREA,Nc,Nr,cW] = define_HSUs(TpInd,TI_bins...
    ,DEMx,B,R,NB,betaD,alphaD,cs,Axs,ch,outletW);
%--------------------------------------------------------------------------
%                              diffusion matrix
%--------------------------------------------------------------------------
if strcmp(DIFFUSION,'on')
  %calculate the diffusion matrix
  D                    = Diffusion_Matrix(HSUs,Nc);
else
  D                    = zeros(Nc,Nc);
end
%--------------------------------------------------------------------------
%                      Flow Distribution Matrices (FDMs)
%--------------------------------------------------------------------------
%multiple-direction DOWNSLOPE FDM for surface
WxmD                 = FDM_calc_surface(MxmD,MxsD,HSUs,Nc,Nr);
%multiple-direction UPSLOPE FDM for surface
if strcmp(DIFFUSION,'on')
  WxmU                 = FDM_calc_surface(MxmU,MxsU,HSUs,Nc,Nr);
else
  WxmU                 = zeros(Nc,Nc);
end
%multiple-direction DOWNSLOPE FDM for subsurface
WbmD                 = FDM_calc_subsurface(MbmD,HSUs,Nc);
%multiple-direction UPSLOPE FDM for subsurface
if strcmp(DIFFUSION,'on')
  WbmU                 = FDM_calc_subsurface(MbmU,HSUs,Nc);
else
  WbmU                 = zeros(Nc,Nc);
end
%--------------------------------------------------------------------------
%visualise
CM=colormap(jet);
CM(1,:)=[];
gry=[1 1 1];

figure(431)
clf
imageschs(DEMx,DEMx.Z,'colormap',CM,'colorbarlabel','Elev','ticklabels','nice','medfilt',1,'nancolor',gry)
drawnow

figure(432)
clf
imageschs(DEMx,TpInd,'colormap',CM,'colorbarlabel','TI','ticklabels','nice','medfilt',1,'nancolor',[1 1 1])
drawnow


figure(433)
clf
imageschs(DEMx,R,'colormap',CM,'colorbarlabel','Reach ID','ticklabels','nice','medfilt',1,'nancolor',[1 1 1])



figure(434)
clf
imageschs(DEMx,HSUs,'colormap',CM,'colorbarlabel','HUSs','ticklabels','nice','medfilt',1,'nancolor',gry)
drawnow


figure(435)
clf
imageschs(DEMx,betaD,'colormap',CM,'colorbarlabel','\beta','ticklabels','nice','medfilt',1,'nancolor',gry)
drawnow

figure(436)
clf
imageschs(DEMx,B,'colormap',CM,'colorbarlabel','isoBs','ticklabels','nice','medfilt',1,'nancolor',gry)
drawnow

%**************************************************************************
function DEMi=Invert_DEM(DEM)
DEMi   = DEM;
dem    = DEMi.Z;
dem    = abs(dem-max(dem(:)));
DEMi.Z = dem;
DEMi   = fillsinks(DEMi);
%**************************************************************************
function D=Diffusion_Matrix(HSUs,Nc)
Edges       = (0.5:1:Nc+0.5)';
D           = zeros(Nc,Nc);
for ii      = 1:Nc
    ix      = HSUs==ii;
    [~,icd] = ixneighbors(HSUs,ix,8);
    A       = HSUs(icd);
    counts  = histcounts(A(:), Edges);     
    D(:,ii) = counts./sum(counts+eps);
end
%**************************************************************************
function [HSUs,SINa,SINb,COSa,COSb,areaf,AREA,Nc,NR,cW]=define_HSUs(TpInd,edgB,DEM,B,R,NB,beta,alpha,cs,As,ch,outletW)
%assign a new ID for reaches in different isobasins to avoid having non unique HSU IDs
PRIMES         = primes(1e8);
PRIMES(1:100)  = [];
Nr             = max(R(:));
for iR         = 1:Nr
  R(R==iR)     = PRIMES(iR);
end
PRIMES(1:Nr+1) = [];

Btemp              = B;
for iB             = 1:NB
  Btemp(Btemp==iB) = PRIMES(iB);
end
%assign each basin a seperate channel segment
R        = R.*Btemp;
%define outlet
[~,ixx]       = sort(As(:),'descend');
outlet        = ixx(1:1);
%exclude outlet from reach map
R(outlet)     = NaN;
%reshuffle to be from 1:NR
R             = shufflelabel(R,true);
%in case NR has changed
NR            = max(R(:));
%calssify hillslope cells
% [~,TPIND0,nc] = classify2(nc,TpInd);
[~,HSUs0,nc] = classify(edgB,TpInd);
%exclude channels: channel reaches will be defined as separate HSUs
HSUs0(ch==1)         = NaN;
%initialise
HSUs         = 0*DEM.Z;
ID            = 0;
%now classify based on TPIND0 and basin ID
for ii = 1:nc
    for jj = 1:NB
        cond   = HSUs0==ii & B==jj;
        
        %ensure condition is met
        if nnz(cond)>0
            ID = ID+1;
            HSUs(cond)=ID;
            %define class width for manning's scaling
            %hilslope classes don't get scaled for width(channel will)
            %their width is equal to DEM resolution
            cW(ID,1)=cs;
        end
    end
end
% add the channel reaches at the end: max(TPIND) will be the outlet reach
%only if we have channels
if NR>0
  %-------------------
  areaMAX=max(As(:));
  for kk = NR:-1:1
    cond = R==kk;
    %ensure condition is met
    if nnz(cond)>0
      ID = ID+1;
      HSUs(cond)=ID;
      %based on Width~Q^0.5
      % we scale them based on outlet: if outlet is 5m then everything else
      % must slightly narrower than the outlet
      cW(ID,1)=outletW*( mean(As(cond)./areaMAX) ).^0.5;
    end
  end
  %-------------------
end
%total number of units
Nc            = max(HSUs(:));
%add one for outlet
Nc            = Nc+1;
%add one for outlet
NR            = NR+1;
%set outlet width
cW(Nc,1)      = outletW;
%update TPIND map
HSUs(outlet) = Nc;
%calculate mean slopes for each class
for ii=1:Nc
    cond          = HSUs==ii;
    b0            = beta(cond);
    a0            = alpha(cond);
    b0(isnan(b0)) = [];
    a0(isnan(a0)) = [];
    SINa(ii,1)    = mean(sind(a0));
    COSa(ii,1)    = mean(cosd(a0));
    SINb(ii,1)    = mean(sind(b0));
    COSb(ii,1)    = mean(cosd(b0));
    areaf(ii,1)   = cs.^2*nnz(cond);
end
%total area
AREA     = sum(areaf);
%fraction of the study area in each TI class
areaf    = (areaf)./AREA;
%impose a 90 degree angle for the outlet cell (because it's a sink)
SINa(Nc) = sind(90);
COSa(Nc) = cosd(90);
SINb(Nc) = sind(90);
COSb(Nc) = cosd(90);

%**************************************************************************
function W = FDM_calc_subsurface(M,TPIND,Nc)
%squeeze HSUs into a vector
TWI     = TPIND(:);
%find the non-zero elements of the flow matrix
[ix,jy] = find(M);
%flow weighting matrix for (HSUs interactions with one another)
W       = zeros(Nc,Nc);
%subsurface isn't bound by the same topographic focusing as the surface
%so it can give flow to other HSUs whether inside the same iso-basin or not
for ii=1:length(ix)
  W( TWI(ix(ii)),TWI(jy(ii)) ) = W( TWI(ix(ii)),TWI(jy(ii)) ) +  M( ix(ii),jy(ii) );
end
%sum of weights should add up to one
SUM        = sum(W,2);
SUM        = repmat(SUM,1,Nc);
W          = W./(SUM+eps);
%last class (catchment outlet) doesn't give to any other HSU
W(end,:)   = 0;
W(end,end) = 1;
W          = W';    
%**************************************************************************
function W = FDM_calc_surface(Mm,Ms,TPIND,Nc,Nr)
%squeeze HSUs into a vector
TWI     = TPIND(:);
%find the non-zero elements of the flow matrix
[ix,jy] = find(Mm);
%flow weighting matrix for (HSUs interactions with one another)
W       = zeros(Nc,Nc);
for ii=1:length(ix)
  %if the giver is a surface hillslope cell
  if TWI(ix(ii))<=Nc-Nr
    %give according to flow MULTIPLE direction distribution matrix
      W( TWI(ix(ii)),TWI(jy(ii)) )=W( TWI(ix(ii)),TWI(jy(ii)) ) +  Mm( ix(ii),jy(ii) );
  else
    %otherwise, if giver is a surface channel cell, then it can only give to
    % other surface channel cells according to SINGLE direction flow matrix 
    if TWI(jy(ii))>Nc-Nr
      W( TWI(ix(ii)),TWI(jy(ii)) ) = W( TWI(ix(ii)),TWI(jy(ii)) ) +  Ms( ix(ii),jy(ii) );
    end
  end
end
%sum of weights should add up to one (for mass continuity)
SUM        = sum(W,2);
SUM        = repmat(SUM,1,Nc);
W          = W./(SUM+eps);
%last class (catchment outlet) doesn't give to any other HSU
W(end,:)   = 0;
W(end,end) = 1;
%transpose for matrix operation within the GMD_TOPMODEL_ode function
W          = W';    
%**************************************************************************
function [Ic,Ind,nc0]=classify(edgB,I)
%classifies based on fixed and uniform bin edges
%mid point values
edg = (edgB(1:end-1)+edgB(2:end))/2;
%loop through
Ic  = 0*I;
Ind = 0*I;
nc0 = 0;
for ii            = 1:length(edg)
    cond          = find(I>=edgB(ii) & I<edgB(ii+1));
    if not(isempty(cond))
        nc0       = nc0+1;
        Ic(cond)  = edg(ii);
        Ind(cond) = nc0;
    end
end
%**************************************************************************
function [Ic,Ind,nc]=classify2(nc,I)
%classifies based on ranks
%sort desxcending
[Is,ix]  = sort(I(:),'descend');
%find NaNs
iNaN     = find(isnan(Is));
%remove NaNs
Is(iNaN) = [];
ix(iNaN) = [];
%define number of cells in each rank
N        = floor(numel(Is)./nc);
%initialise
Ic       = 0*I;
Ind      = 0*I;
%loop through number classes
ID       = nc+1;
for ii            = 1:nc
    if ii<nc
        cond      = (ii-1)*N+1:ii*N;
    else
        cond      = (ii-1)*N+1:numel(Is);
    end
    ID            = ID-1;
    Ic(ix(cond))  = mean(I(ix(cond)));
    Ind(ix(cond)) = ID;
end
