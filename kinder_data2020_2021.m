function kinder_data2020_2021

path      = 'C:\Users\Sal\Desktop\USB\Salim_Newcastle Uni\Reports\Paper 4-catchment scale\study2\DATA\';


% data=xlsread([path 'Kinder_2020_21.xlsx'],'firm_');

data=xlsread([path 'Kinder_master_2020_21.xlsx'],'nog_');


% yyyymmddHH0=[2020 8 22 0];

obsQ=data(:,2)./1000;
obsR=data(:,3)/1000;
obsQ(isnan(obsQ))=eps;
obsR(isnan(obsR))=eps;
DT=(600:600:length(obsQ)*600)';
DTR=DT;



DT0=DT/60/60/24;


cond=DT0<170 | DT0>230;
obsQ(cond)=[];
obsR(cond)=[];
DT=(600:600:length(obsQ)*600)';
DTR=DT;
yyyymmddHH0=[2021 06 19 0];



save([path 'obsData' '2021' 'nog_' 'calib'],'DT','obsQ','DTR','obsR','yyyymmddHH0');




