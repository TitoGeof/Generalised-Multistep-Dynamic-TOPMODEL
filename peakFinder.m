function [oQpeaks,oQtimes,pQpeaks,pQtimes,Qfrac] = peakFinder(NPeaks,oQ,pQ,QFRAC,dTime)
%left/right window to look for peaks 
C       = 10/(dTime/60);
%find peaks in observed discharge
[oQpeaks,oQtimes,~,promO] = findpeaks(oQ,'MinPeakDistance',8*C);
%sort according to prominence
[~,ixo] = sort(promO,'descend');
%choose the top NPeaks peaks
ixo     = ixo(1:NPeaks);
oQpeaks = oQpeaks(ixo);
oQtimes = oQtimes(ixo);
%define windows before and after the obs data peak locations
window  = 12*C;
pQtimes = [];
Qfrac   = [];
for ii  = 1:NPeaks
    %define window index
    tleft  = oQtimes(ii)-24*C;   tleft(tleft<1)=1;
    tright = oQtimes(ii)+24*C;   tright(tright>length(pQ)) = length(pQ);
    %find the peaks in predicted discharge
    [peaks,pQt_id,~,prom] = findpeaks(pQ(tleft:tright));
    %check if enough peaks are identified
    if not(isempty(peaks))
        %ensure to not pick one peak two times
        for kk = 1:length(pQtimes)
            id         =find(pQt_id+tleft-1==pQtimes(kk));
            peaks(id)  = [];
            pQt_id(id) = [];
            prom(id)   = [];
        end
        if not(isempty(peaks))
            %sort according to prominence
            [~,ixp]       = sort(prom,'descend');
            ixp           = ixp(1);
            pQpeaks(ii,1) = peaks(ixp);
            pQtimes(ii,1) = pQt_id(ixp)+tleft-1;
            Qfrac(ii,1)   = QFRAC(pQt_id(ixp)+tleft-1);
        else
            pQpeaks(ii,1) = 0;
            pQtimes(ii,1) = -1;
            Qfrac(ii,1)   = 1;
        end
    else
        pQpeaks(ii,1) = 0;
        pQtimes(ii,1) = -1;
        Qfrac(ii,1)   = 1;
    end
end