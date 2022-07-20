function [obsQ,obsR,DATETIME,yymmddHH0]=LoadObsData(year,catchname,period)
% DATETIME=readtable([year '_refinedData_meanRain.xlsx'],'sheet',catchname);
DATETIME=readtable([year '_refinedData.xlsx'],'sheet',catchname);
DATETIME(:,2:end)=[];
%convert table format to array
DATETIME=table2array(DATETIME);
if iscell(DATETIME)==1
    DATETIME=string(DATETIME);
    if strlength(DATETIME(2))>17
        DATETIME = datetime(DATETIME,'InputFormat','dd/MM/yyyy HH:mm:ss');
    else
        DATETIME = datetime(DATETIME,'InputFormat','dd/MM/yy HH:mm:ss');
    end
end

%for Nosong 2010, shoft the datetime one hour
if strcmp(catchname,'nog_')==1 && strcmp(year,'2010')==1
    shift=DATETIME(7)-DATETIME(1);
    DATETIME=DATETIME-shift;
end

% data=xlsread([year '_refinedData_meanRain.xlsx'],catchname);
data=xlsread([year '_refinedData.xlsx'],catchname);

if strcmp(year,'2010')==1
    if strcmp(period,'calib')==1
        %calibration period
        yymmddHH0=[2010 08 22 00];
        DT0=datetime([2010 08 22 00 00 00]);
        DTN=datetime([2010 10 10 00 00 00]);
    elseif strcmp(period,'valid')==1
        %validation period
        yymmddHH0=[2010 10 08 00];
        DT0=datetime([2010 10 08 00 00 00]);
        DTN=datetime([2010 11 15 00 00 00]);
    end
elseif strcmp(year,'2012')==1
    if strcmp(period,'calib')==1
        %calibration period
        yymmddHH0=[2012 08 22 00];
        DT0=datetime([2012 08 22 00 00 00]);
        DTN=datetime([2012 10 10 00 00 00]);
    elseif strcmp(period,'valid')==1
        %validation period
        yymmddHH0=[2012 10 08 00];
        DT0=datetime([2012 10 08 00 00 00]);
        DTN=datetime([2012 11 15 00 00 00]);
    end
end
ind=DATETIME<=DTN & DATETIME>=DT0;
obsQ=data(ind,2);
obsR=data(ind,1);
DATETIME=DATETIME(ind);
%convert Lit/s to m3/s
obsQ=obsQ*1e-3;
%convert mm/s to m/s
obsR=obsR*1e-3;