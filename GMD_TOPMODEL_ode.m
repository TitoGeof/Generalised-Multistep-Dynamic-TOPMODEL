function [Qt,Qfrac,simTime,dTime,massErr,KGE,SpinUp]=GMD_TOPMODEL_ode(D,WxmD,WxmU,WbmD,WbmU,obsR,obsQ,...
                        params,cs,DT0,t,SINa,SINb,COSa,COSb,area,AREA,Nr,Nc,cW,SpinUp)              
warning('off','all');
%-------------------------------------------------------------------------- 
%number of observed data points [-]
Nobs                 = size(obsR,1);
%observed rainfall resolution [s]
dTime                = t(2)-t(1);
%define spin up in terms of data timestep
SpinUp               = SpinUp*(24*60*60/dTime);
%load uncertain/input model parameters
[d,Tmax,ep,Smax,mannNhs,mannNch,Hmax,ABSTOL] = unPack_uncertain_parameters(params);
%initialise system
[Sx0,Su0,Sw0]        = initialiseSYS(Nc,Hmax);
%assemble manning's n coefficient for hillslope vs channel classes
mannN                = 0*Sx0 + mannNhs;           %hillslope
mannN(Nc-Nr+1:end,1) = mannNch;                   %channels
%annual daily average evapo-transpiration rate [m/s]
ET                   = seasonal_sinewave_evap(DT0,1,dTime,(1:length(obsR))');
%convert total rainfal (e.g., from tipping bucket) to rainfall intensity [m/s]
drdt                 = obsR./dTime;
%define method of interpolation
METHOD               = 'pchip';
%create gridded interpolation
FDR                  = griddedInterpolant(t,drdt,METHOD);
FET                  = griddedInterpolant(t,ET  ,METHOD);
%assemple the vector of state variables
V0                   = [Sx0;Su0;Sw0];
%--------------------------------------------------------------------------
%jacobian pattern (tells the ode-solver to evaluate only where JPAT=1 to save runtime - see model paper for details)
M1                   = full(diag(ones(Nc,1)));
M2                   = sign(WxmD + D' + WxmU);
M3                   = sign(WbmD + D' + WbmU);
M2(M1==1)            = 1;
M3(M1==1)            = 1;
JPAT                 = [M2 M1 M1; M1 M1 M1; M1 M1 M3];      
%--------------------------------------------------------------------------
%define ode-solver time vector
tO                   = linspace(t(1),t(end),5*length(t))';
%ode-solver Options
OPS                  = odeset('JPattern',JPAT,'InitialStep',1e-64,'AbsTol',ABSTOL,'RelTol',100*ABSTOL);
%solve using ode15s, suitable for "stiff" system if equations
tic;
[~,V]                = ode15s(@HydroGEM_ode_fun,tO,V0,OPS,area,d,Nc,Smax,FET,FDR,mannN,D,WxmD...
                              ,WbmD,WxmU,WbmU,cs,SINa,SINb,COSa,COSb,Tmax,Hmax,Nr,cW,ep);
simTime              = toc;
e                    = 1e-64;
%in case ode-solver has crashed midway
if size(V,1)<Nobs; V = nan(Nobs,3*Nc); end
%--------------------------------------------------------------------------
%calculate outlet discharge
%--------------------------------------------------------------------------
%disaggregate variables 
Sx                   = V(:,Nc);
Sw                   = V(:,3*Nc);
%calculate base flow
T                    = Tmax.*(Sw./Hmax).^d;
qb                   = SINa(Nc).*T./cs;
qb(qb>Sw)            = Sw(qb>Sw);
qb                   = area(Nc).*qb*AREA;
%remove rootzone storage from total Sx to caluculate overland flow 
Sx                   = Sx-Smax;
Sx(Sx<0)             = 0;
%hydraulic radius for channel class
Ss                   = Sx./cW(Nc)*cs;
R                    = Ss.*cW(Nc)./(2*Ss+cW(Nc));
%manning's velocity
v                    = R.^(2/3).*SINb(Nc).^(1/2)./mannN(Nc);
qo                   = Sx.*v/cs;
qo(qo>Sx)            = Sx(qo>Sx);
qo                   = area(Nc).*qo*AREA;
%base flow fraction
Qfrac                = qb./(qo+qb+e);
%calculate total flow (hydrograph)
Qt                   = qo + qb;
%--------------------------------------------------------------------------
%calculate mass error
[massErr]            = massError(tO,Nc,Smax,ep,V,Qt,FDR,FET,AREA,area,SpinUp);
%--------------------------------------------------------------------------
%interpolate predicted discharge to observed discharge/rainfall times
Qt                   = interp1(tO,Qt,t,'pchip');
Qfrac                = interp1(tO,Qfrac,t,'pchip');
%--------------------------------------------------------------------------
%calculate Kling Gupta performance metric
[KGE]                = ObjectiveFunCalculation(Qt,obsQ,SpinUp);
%**************************************************************************
%**************************************************************************
%                               subfunctions
%**************************************************************************
%**************************************************************************
%**************************************************************************
function dVdt=HydroGEM_ode_fun(t,V,area,d,Nc,Smax,FET,FDR,mannN,D,WxmD...
                              ,WbmD,WxmU,WbmU,cs,SINa,SINb,COSa,COSb,Tmax,Hmax,Nr,cW,ep)    
%--------------------------------------------------------------------------                                                
%machine precision
e                = 1e-64;
%--------------------------------------------------------------------------
%                          disaggragte variables
%--------------------------------------------------------------------------
V                = real(V);
Sx               = V(1:Nc,1)  ;          
Su               = V(Nc+1:2*Nc,1) ;     
Sw               = V(2*Nc+1:3*Nc,1) ; 
%--------------------------------------------------------------------------
%          interpolate rainfall and evapo-transpiration rates
%--------------------------------------------------------------------------
%interpolate rainfall intensity based on solver time
Rn               = FDR(t);   Rn(Rn<e) = e;
%interpolate potential evapotranspiration rate
ET               = FET(t);  
Ep               = ep.*ET;   Ep(Ep<e)=e;
%--------------------------------------------------------------------------
%                       subsurface storage updates
%--------------------------------------------------------------------------
%fraction from saturated zone to the surface
w0               = stepfun(Sw,Hmax,e);
qwt              = w0.*(Sw-Hmax);
Sw               = Sw - qwt;
%fraction from unsaturated zone to the surface
Hu               = Hmax-Sw;
w1               = stepfun(Su,Hu,e);
quz              = w1.*(Su-Hu);
Su               = Su - quz;
%--------------------------------------------------------------------------
%                    water table inflows and outflows
%--------------------------------------------------------------------------
%subsurface power-law transmissivity profile
T                = Tmax.*(Sw./Hmax).^d;
%subsurface diffusion
dSwi             = repmat(Sw',Nc,1)-Sw;
dSwdx            = sum(D.*dSwi,2)/cs;
%diffusion flux
flxb             = SINa - COSa.*dSwdx;
%determine the sign (+:downslope, -:upslope)
wb               = stepfun(flxb,0,e);
%downslope flow velocity
qbD              = wb.*flxb.*T./cs;
%upslope flow veolcity
qbU              = (1-wb).*abs(flxb).*T./cs;
%route downslope and upslope flows
[qbiD,qbD]       = distribute_flow(WbmD,area,Sw,qbD,Nc,e); 
[qbiU,qbU]       = distribute_flow(WbmU,area,Sw,qbU,Nc,e);
%combine downslope and upslope
qb               = qbD  + qbU;
qbi              = qbiD + qbiU;
%--------------------------------------------------------------------------
%                  unsaturated zone inflows and outflows
%--------------------------------------------------------------------------
%vertical hydraulic conductivity averaged across unsaturated zone's thickness
Kbar              = Su.*(Tmax-T)./(Hu+e).^2;
%total vertical flow based on Richard's type Equation
qv                = Kbar.*(-Su./(Hu+e).*COSb+1);  
%ensure it doesn't exceed available Su
[w2,qv]           = stepfun(qv,Su,e);
qv                = w2.*Su+(1-w2).*qv;
%--------------------------------------------------------------------------
%                         surface storage updates
%--------------------------------------------------------------------------
%calculate actual evapo-transpiration
w3                = stepfun(Sx,Smax,e);
Srz               = w3.*Smax + (1-w3).*Sx;
Ea                = Ep.*Srz./(Smax);
[w4,Ea]           = stepfun(Ea,Srz,e);
Ea                = w4.*Srz + (1-w4).*Ea;
%update surface excess storage available for routing
Sx                = Sx - Srz ;   
%total available storage in the subsurface (saturated + unsaturated)
SD                = Hmax-Su-Sw;
%in case of no infiltration excess, give all you can to the subsurface
w5                = stepfun(Sx,SD,e);
qx                = w5.*SD + (1-w5).*Sx;
%--------------------------------------------------------------------------
%                       surface inflows and outflows
%--------------------------------------------------------------------------
%scale surface excess to acount for variable channel width (only in channel HSUs)
%hilslope cW is set equal to cs, so no scaling occurs there
Ss                = Sx*cs./cW;
%hydraulic radius for channel HSUs (assuming rectangular channel)
R                 = Ss.*cW./(2*Ss+cW);
%hydraulic radius for hillslope HSUs
R(1:Nc-Nr)        = Ss(1:Nc-Nr);
%overland flow diffusion
dSxi              = repmat(Sx',Nc,1)-Sx;
dSxdx             = sum(D.*dSxi,2)./cs;
%diffusion flux
flxx              = SINb - COSb.*dSxdx;
%determine the sign (+:downslope, -:upslope)
wx                = stepfun(flxx,0,e);
%Manning's n downslope flow flux (v*depth/grid resolution)
qoD               = wx.*R.^(2/3).*flxx.^(1/2)./mannN.*Sx./cs;
%upslope flow flux 
qoU               = (1-wx).*R.^(2/3).*abs(flxx).^(1/2)./mannN.*Sx./cs;
%route downslope and upslope flows
[qoiD,qoD]        = distribute_flow(WxmD,area,Sx,qoD,Nc,e);
[qoiU,qoU]        = distribute_flow(WxmU,area,Sx,qoU,Nc,e);
%combine downslope and upslope
qo                = qoD  + qoU;
qoi               = qoiD + qoiU;
%--------------------------------------------------------------------------
%                       update all time derivatives 
%--------------------------------------------------------------------------
dSxdt             = qoi - qo - qx + quz + qwt + Rn - Ea; %surface storage
dSudt             = qx  - qv - quz;                      %unsaturated zone 
dSwdt             = qbi - qb + qv - qwt;                 %subsurface storage
%aggregate derivatives for the ODE-solver
dVdt              = [dSxdt;dSudt;dSwdt];
%**************************************************************************
function [qoi,qo] = distribute_flow(W,area,S,q,Nc,e)
%ensure it does not exceed the available storage
[w7,q]            = stepfun(q,S,e);
qo                = (1-w7).*q + w7.*S;
%flux into units
qq                = area.*qo;
qq(Nc)            = 0;% last HSU (outlet) does not give flow to other HSUs
qoi               = (W*qq)./(area+eps);
%**************************************************************************
function [F,x]    = stepfun(x,x0,e)
%a step function to handle the discontinous solution profiles
x(x<e)            = e;
F                 = (1+tanh((x-x0)./e))/2;
%**************************************************************************
function [KGE]=ObjectiveFunCalculation(pQ,oQ,SpinUp)
%exclude Spin-up period from evaluaiton
oQ(1:SpinUp)     = [];
pQ(1:SpinUp)     = [];
%--------------------------------------------------------------------------
%KlingGupta Efficiency
KGE = 1-sqrt( (corr(pQ,oQ)-1).^2 + (std(pQ)./std(oQ)-1).^2 + (mean(pQ)./mean(oQ)-1).^2 );
%**************************************************************************
function [massErr]=massError(t,Nc,Smax,ep,V,pQ,FDR,FET,AREA,area,SpinUp)
%machine precision
e                = 1e-64;
%exclude Spin-up period from evaluaiton
t(1:SpinUp)      = [];
V(1:SpinUp,:)    = [];
pQ(1:SpinUp)     = [];
%--------------------------------------------------------------------------
%data timestep 
dTime=t(2)-t(1);
%--------------------------------------------------------------------------
%interpolate rainfall intensity based on solver time
Rn               = FDR(t);   Rn(Rn<e) = e;
%calculate volume of rain [m3]
RnT              = Rn*AREA*dTime;
%--------------------------------------------------------------------------
%interpolate potential evapotranspiration rate
ET               = FET(t);   ET(ET<e) = e;
Ep               = ep.*ET;
%amount of surface water
Sx               = V(:,1:Nc);
%calulate rootzone storage available for evapo-transpiration
w3                = stepfun(Sx,Smax,e);
Srz               = w3.*Smax + (1-w3).*Sx;
%calculate actual evapo-transpiration
Ea                = Ep.*Srz./(Smax);
[w4,Ea]           = stepfun(Ea,Srz,e);
Ea                = w4.*Srz + (1-w4).*Ea;
%calculate total Ea across the catchment in each timestep
EaT               = sum(area'.*Ea*AREA*dTime,2);
%--------------------------------------------------------------------------
%calculate total storage across the catchment in each timestep
A                 = [area;area;area]'*AREA;
ST                = sum(A.*V,2);
%--------------------------------------------------------------------------
%mass balance in each timestep as a percentage of total rain
%*** |dS - r + Ea + Q| = 0 in each timestep ***
err              = abs(ST-ST(1)- cumsum(RnT) + cumsum(EaT) + cumsum(pQ*dTime));
ERR              = err./(sum(RnT(:))+1)*100;
%error at the end of simulation
massErr          = ERR(end);
