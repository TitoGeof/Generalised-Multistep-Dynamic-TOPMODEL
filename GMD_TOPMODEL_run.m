function GMD_TOPMODEL_run
clc;
%--------------------------------------------------------------------------
%Generalised Multistep Dynamic(GMD) TOPMODEL is described in detail in a paper 
%with the same title, published in Water Resources Research Journal, in 2022.
%Developer: Dr Salim Goudarzi: salim_goodarzi@yahoo.com
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
catchname   = 'nog12_'; %cacthment file name (your DEM should be a .TIFF file called [catchname 'dem.tiff'])
%and all data files should be placed in a folder called DATA in the
%directory where the model files are, i.e., data files to be loaded by the 
% model should be in path=[pwd '\DATA\'].
%--------------------------------------------------------------------------
DIFFUSION   = 'off'; % 'on' or 'off'. DIFFUSION: if 'on' diffusion flow routing
%is performed on both hillslope and channel HSUs. if 'off', kinematic wave
%routing is excecuted (if turned on or off, the "Preprocess_DEM" routine needds
%to be performed again, to make sure diffusion matrix is included/excluded)
%--------------------------------------------------------------------------
Href        = 0.5; %reference head-drop along flow path path [m] (for calculating 
%subsurface hydraulic gradient) it is catchment-specific. Refer to the model
%jouranl article for more info.
%--------------------------------------------------------------------------
%number Topographic Index (TI) classes
nClass       = 7;
%--------------------------------------------------------------------------
ISOBASINS    = 'on';   % 'on' or 'off'. if 'off' discretisation is done based only on
%Topographic Index (TI), if 'on'`, discretisation is done by the combination of TI
%and iso-basins of area ~ 'Athresh'
%--------------------------------------------------------------------------
Athresh      = 2000; %iso-basin area threshold [m^2]: this value should be
%set according to your application but also catchment size. The iso-basin
%code can take a long time to run if this number is small, or catchment is
%very large; worse if both. For very large catchments (>100km2) it may take
%a few days. But it is part of the DEM pre-processing and only needs to run once.
%--------------------------------------------------------------------------
CHthresh     = 200; %channel initiation area threshold [m^2]
outletW      = 2;   %catchment outlet channel width [m]
%--------------------------------------------------------------------------
cs           = 2; %DEM cellsize [m]; resolution of DEM raster
%--------------------------------------------------------------------------
SpinUp       = 10; %a 'warm-up' period (predictions discarded here)  
%--------------------------------------------------------------------------
%                 input model parameters (need calibration)
%--------------------------------------------------------------------------
%power law exponent of conductivity decay with depth, d, [m]
PARAMset(1)  = 8;
%maximum transmissivity (at saturation), Tmax [m2/s]
PARAMset(2)  = 3e-4;
%effective frainable porosity, phi, [m]
PARAMset(3)  = 4e-3;
%maximum daily evaporation rate, averaged across a year, ep [m/day]
PARAMset(4)  = 0.0037;
%maximum root-zone storage, Smax [m]
PARAMset(5)  = 0.01;
%hilslope Manning's n_{hs} [s/m^(1/3)] 
PARAMset(6)  = 1;
%channels Manning's n_{ch} [s/m^(1/3)]
PARAMset(7)  = 0.5;
%--------------------------------------------------------------------------
%                  read catchment rainfall and discharg data
%--------------------------------------------------------------------------
load([pwd '\DATA\' 'obsData' catchname],'DT','obsQ','DTR','obsR','yyyymmddHH0');
%obsQ [m3/s]: observe discharge record
%obsR [m]: obsrved rainfall record
%DT: time of observed discharge. It has to be in seconds.
%DTR: time of observed rainfall. It has to be in seconds.

%NOTE: that obsR and obsQ don't need to have the same resolution necessarily
%because GMD-TOPMODEL is time-continuous model. However, make sure obsR and obsQ
%have the correct DT and DTR vectors if they are at different time scales.

%yyyymmddHH0: is a four element array, marking the staring date-time of your
%simulation (record). It is used to calculate seasonal fluctuations in 'ep' 
%yyyymmddHH0(1):year, yyyymmddHH0(2):month, %yyyymmddHH0(3):day, yyyymmddHH0(4):hour
%--------------------------------------------------------------------------
%                       pre-processing DEM data
%--------------------------------------------------------------------------
%NOTE: run this function only once and save outputs, because it can take a 
%long time for large catchments. Especially if ISOBASINS and/or DIFFUSION
%are turned 'on'. After that, load values from disk for future simulations of you catchment.

[WxmD,WxmU,WbmD,WbmU,D,SINa,SINb,COSa,COSb,areaf,AREA,Nc,Nr,TPIND,cs,DEM,cW]...
     =Preprocess_DEM(catchname,CHthresh,Href,Athresh,cs,ISOBASINS,DIFFUSION,nClass,outletW,PARAMset(1));

NAME = [catchname num2str(Href) '_' num2str(Athresh) 'DEM_PRE_PROCESS' '_d=' num2str(PARAMset(1))];
save([pwd '\DATA\' NAME '.mat'],'areaf','AREA','WxmD','WxmU','WbmD','WbmU','SINa','SINb','COSa','COSb','cs','TPIND','Nr','Nc','DEM','D','cW')
load([pwd '\DATA\' NAME '.mat'],'areaf','AREA','WxmD','WxmU','WbmD','WbmU','SINa','SINb','COSa','COSb','cs','TPIND','Nr','Nc','DEM','D','cW')
%--------------------------------------------------------------------------
%                              run GMD-TOPMODEL 
%--------------------------------------------------------------------------
[predQ,Qfrac,simTime,dTime,KGE,SpinUp] = GMD_TOPMODEL_ode(D,WxmD,WxmU,WbmD,WbmU,obsR,obsQ...
  , PARAMset,cs,yyyymmddHH0,DTR,SINa,SINb,COSa,COSb,areaf,AREA,Nr,Nc,cW,SpinUp);
%--------------------------------------------------------------------------
%                         plot rainfall-discharge 
%--------------------------------------------------------------------------
DT      = DT/60/60/24;
DTR     = DTR/60/60/24;
obsR    = obsR./dTime*60*60*1000;


TITLE=['catchment: ' num2str(round(AREA/1e6,1,'significant')) 'km^2' ' | ' '#HSUs: ' num2str(Nc) ' | ' 'runtime: ' num2str(round(simTime)) 's' ...
  ' | ' 'KGE: ' num2str(round(KGE*100)) '%'];
%--------------------
figure(201)
clf
%--------------------
hh1=subplot(3,1,1:2);
%--------------------
title(TITLE)
hold on
[haxes1,hline1,hline2] = plotyy(DT,obsQ,DTR,obsR,'area','area');
plot(DT,predQ,'k:','linewidth',1.5)
plot(DT,Qfrac.*predQ,'r-')
plot([DT(SpinUp) DT(SpinUp)],[0 max(obsQ)*100],'m--','linewidth',1)
ylabel(haxes1(1),'Q [m^3/s]','color','k');
set(hline1(1),'FaceColor','c','EdgeColor','c');
set(haxes1(1),'YColor','k')
ylim(haxes1(1),[min(Qfrac.*predQ) max(obsQ)*1.5]);
xlim(haxes1(1),[min(DT) max(DT)]);
ylabel(haxes1(2),'Rain [mm/hr]','color','b');
set(hline2,'FaceColor','b','EdgeColor','b');
set(haxes1(2),'YColor','b');
set(haxes1(2),'YDir','reverse');
ylim(haxes1(2),[0 max(obsR)*5]);
xlim(haxes1(2),[min(DT) max(DT)]);
LEG=legend('obsQ','predQ','baseflow'...
,'spin-up','rainfall','location','best','orientation','vertical');
set(LEG,'box','on','color','w')
%--------------------
hh2=subplot(313);
%--------------------
hold on
[haxes1,hline1,hline2] = plotyy(DT,obsQ,DTR,obsR,'area','area');
plot(DT,predQ,'k:','linewidth',1.5)
plot(DT,Qfrac.*predQ,'r-')
plot([DT(SpinUp) DT(SpinUp)],[min(Qfrac.*predQ) max(obsQ)*100],'m--','linewidth',1)
ylabel(haxes1(1),'Q [m^3/s]','color','k');
ylim(haxes1(1),[1e-6 max(obsQ)*1.5]);
xlim(haxes1(1),[min(DT) max(DT)]);
set(hline1(1),'FaceColor','c','EdgeColor','c');
set(haxes1(1),'YColor','k')
ylabel(haxes1(2),'Rain [mm/hr]','color','k');
set(haxes1(2),'YDir','reverse');
ylim(haxes1(2),[0 max(obsR)*10]);
xlim(haxes1(1),[min(DT) max(DT)]);
set(haxes1(2),'YColor','b');
set(hline2,'FaceColor','b','EdgeColor','b');
set(gca,'yscale','log')
xlabel(haxes1(1),'Time [day]','color','k');
linkaxes([hh1 hh2],'x')
