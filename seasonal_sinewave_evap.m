function Ep=seasonal_sinewave_evap(yyyymmddHH,Ep0,dTime,iTime)
%disagregate the start time/date of the simulation
yyyy0        = yyyymmddHH(1);
mm0          = yyyymmddHH(2);
dd0          = yyyymmddHH(3);
HH0          = yyyymmddHH(4);
%define the starting timedate
datetime0    = datetime(yyyy0,mm0,dd0);
%define 1st of january at 00:00:00 of that year as the reference datetime
datetime_ref = datetime(yyyy0,01,01);
%calculate the number of days from the reference datetime
dayNum       = datenum(datetime0)-datenum(datetime_ref) + floor(iTime.*dTime/24/60/60);
%ensure dayNum is not greater than 365
dayNum       = dayNum-floor(dayNum/365)*365;
%calculate the fraction based on the day of the year
dayFrac      = 0.5*(2+sin(2*pi*dayNum/365-3*pi/4));
% %calculate hour of the day
hr           = HH0+floor(iTime*dTime/60/60);
% %ensure the hour value is not greater than 24hr
hr           = hr-floor(hr/24)*24;
hrFrac       = 0.5*(2+sin(2*pi*hr/24-pi/2));
hrFrac       = hrFrac./max(hrFrac);
Ep           = Ep0*dayFrac.*hrFrac;