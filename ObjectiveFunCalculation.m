%**************************************************************************
function [ofsALL,oQtimes,oQpeaks,pQtimes,pQpeaks]=ObjectiveFunCalculation(...
    SpinUp,pQ,oQ,QFRAC,dTime)
%remove spin up period from analysis (normally first few days at the beginning of the record)
oQ(1:SpinUp)    = [];
pQ(1:SpinUp)    = [];
QFRAC(1:SpinUp) = [];
%calculates top N peaks' magnitudes and timings (perfromance metrics)
NPeaks = 3;
[oQpeaks,oQtimes,pQpeaks,pQtimes,Qfrac] = peakFinder(NPeaks,oQ,pQ,QFRAC,dTime);
%mean peak magnitude error (% of max peak)
peakMagnErr = abs(oQpeaks-pQpeaks)./oQpeaks*100;
ofs_RMQ     = mean(peakMagnErr);
%remove cases where peak did not exist (pQtimes==-1)
oQtimes0      = oQtimes;
pQtimes0      = pQtimes;
cond          = pQtimes==-1;
pQtimes(cond) = [];
oQtimes(cond) = [];
Qfrac(cond)   = [];
if isempty(pQtimes)
    pQtimes   = 0;
    oQtimes   = 0;
end
%mean peak timing error (mins)
peakTimeErr   = abs(pQtimes - oQtimes )*dTime /60;
ofs_RMT       = mean(peakTimeErr);
%for plotting
oQtimes       = oQtimes0;
pQtimes       = pQtimes0;
pQtimes(cond) = oQtimes(cond);
%objective functions: Nash-Sutcliffe efficiency (high and low flows)
pQ            = pQ+1e-8;
oQ            = oQ+1e-8;
ofs_NSE = 1-sum( (oQ-pQ).^2 )./sum( (oQ-mean(oQ)+eps).^2 );
ofs_KGE = 1-sqrt( (corr(pQ,oQ)-1).^2 + (std(pQ)./std(oQ)-1).^2 + (mean(pQ)./mean(oQ)-1).^2 );
%aggregate objective functions & performance metrics
ofsALL = [ofs_NSE,ofs_KGE, ofs_RMQ, ofs_RMT];
