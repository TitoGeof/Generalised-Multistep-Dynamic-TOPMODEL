function[W,w]=distribution_around_mean(cond,predQ_all)
predQ_all=predQ_all+0.1;
%epsilon
e=1e-8;
%store all predQs
pQ=predQ_all;
%remove nonbehavs
pQ(cond,:)=[];
%calculate arithmetic mean
pQm=mean(pQ,1)';
%total number of data points in the record
Nobs=length(pQm);
%pre allocate array
W=zeros(size(predQ_all,1),1);
%calculate standard error
w=objfuns(pQ',pQm)';
%weight according to lowest standard error
W(1-cond==1)=w;


function w=objfuns(pQ,pQm)

ofs_NSE_H = 1-sum( (pQm-pQ).^2 )./sum( (pQm-mean(pQm)+eps).^2 );
% ofs_MAE = mean( abs(oQ-pQ) );
ofs_NSE_L = 1-sum( (log(pQm)-log(pQ)).^2 )./sum( (log(pQm)-log(mean(pQm))).^2);
% ofs_VE = 1-sum( abs(oQ-pQ) )./sum( (oQ+eps) );
% ofs_NSE = 1-sum( abs(oQ-pQ) )./sum( abs(oQ-mean(oQ)+eps) );


ofs_RSME_H = sqrt(sum( (pQm-pQ).^2 )./length(pQm));
ofs_RSME_L = (((1+pQ).^0.3-1)./0.3 -((1+pQm).^0.3-1)./0.3).^2;
ofs_RSME_L= sqrt(sum(ofs_RSME_L)./length(pQm));

[w_NSEH]=ofs_weight(ofs_NSE_H,'NSE');
[w_NSEL]=ofs_weight(ofs_NSE_L,'NSE');
[w_RMSEH]=ofs_weight(ofs_RSME_H,'NONE');
[w_RMSEL]=ofs_weight(ofs_RSME_L,'NONE');

w=w_NSEH.*w_NSEL.*w_RMSEH.*w_RMSEL;

w=w./sum(w+eps);


%**************************************************************************
function [w,ofs]=ofs_weight(ofs,CASE)
if strcmp(CASE,'NSE')
    worst=min(ofs);
else
    worst=max(ofs);
end
ofs(isnan(ofs))=worst;
w=sqrt((ofs-worst).^2);
% w=abs(ofs-worst);
w=w./(max(w)+eps);


