function reduce_sets
clc
warning('off','all');

path='C:\Users\Sal\Desktop\USB\Salim_Newcastle Uni\Reports\Paper 2-time continuous\study2\DATA\';

nSims     = 10000;
catchname = 'nog_';
year      = '2010';
deltaT    = {'10min'};
period    = 'calib';
batch     = '(3)';
MODEL     = {'dynamic','dynamic_ode'};

for idt=1:length(deltaT)
    for K=1:length(MODEL)
        %load relevant parameter-set
        load([path 'PARAMset' batch '_'  num2str(nSims) '.mat'],'PARAMset');
        %predicted discharge
        filename=[MODEL{K} year '_' period '_' catchname num2str(nSims) '_' batch deltaT{idt} ];
        load([path 'ofs_' filename '.mat'],'ofsALL');
        %calculate weights for all objective functions
        [~,w,ofsALL,Nb,del]=LikelihoodWeightCal(ofsALL);
        load([path filename '.mat'],'predQ_all','QFRAC');
        predQ_all(del,:)=[];
        QFRAC(del,:)=[];
        PARAMset(del,:)=[];
        if strcmp(MODEL{K},'dynamic')
            save([path filename 'TA@10min'],'predQ_all','QFRAC', '-v7.3');
            save([path 'ofs_' filename 'TA@10min'],'ofsALL','w','Nb','-v7.3');
            save([path 'PARAMset' batch '_'  num2str(nSims) MODEL{K} 'TA@10min' '.mat'],'PARAMset')
        else
            save([path filename 'TC@10min'],'predQ_all','QFRAC', '-v7.3');
            save([path 'ofs_' filename 'TC@10min'],'ofsALL','w','Nb','-v7.3');
            save([path 'PARAMset' batch '_'  num2str(nSims) MODEL{K} 'TC@10min' '.mat'],'PARAMset')
        end
        
        disp(filename)
    end
end

