%**************************************************************************
function [W,w,ofsALL,Nb,cond]=LikelihoodWeightCal(ofsALL)
%--------------------------------------------------------------------------
%                   evaluate weights for all simulations
%--------------------------------------------------------------------------
W=zeros(length(ofsALL(:,1)),1);
%remove sets outside LOA and OFL
cond=ofsALL(:,4)==0 | ofsALL(:,5)==0 ;
ofsALL(cond,:)=[];
%number of behavioural simulations
Nb=size(ofsALL,1);
%unpack
ofs_NSEH=ofsALL(:,1);
ofs_RMQ=ofsALL(:,2) ;
ofs_RMT=ofsALL(:,3) ;
ofs_QFR=ofsALL(:,4) ;
ofs_LOA=ofsALL(:,5) ;
ofs_NSEL=ofsALL(:,6);
ofs_KGEH=ofsALL(:,7);
ofs_KGEL=ofsALL(:,8);
%calculate performacen metric weights
[w_NSEH]=ofs_weight(ofs_NSEH,'NSE');
[w_RMQ]=ofs_weight(ofs_RMQ,'RMQ');
[w_RMT]=ofs_weight(ofs_RMT,'RMT');
[w_NSEL]=ofs_weight(ofs_NSEL,'NSE');
[w_KGEH]=ofs_weight(ofs_KGEH,'KGE');
[w_KGEL]=ofs_weight(ofs_KGEL,'KGE');


%decide which to include in weighting process by multiplying
% w=w_NSE.*w_RMQ.*w_RMT ;
% w=w_NSE;
% w=w_NSEH.*w_NSEL.*w_RMQ.*w_RMT;
w=w_NSEH.*w_NSEL;
% w=ofs_NSEH;
% w=w_NSEH.*w_NSEL.*w_KGEH.*w_KGEL;


% [ww,ix]=sort(w,'descend');
% w(ix(50:end))=0;

w=w./sum(w+eps);
%replace the weights in the right place
W(1-cond==1)=w;
% W=w;
%**************************************************************************
function [w,ofs]=ofs_weight(ofs,CASE)
if strcmp(CASE,'NSE') || strcmp(CASE,'KGE') 
    worst=min(ofs);
else
    worst=max(ofs);
end
ofs(isnan(ofs))=worst;
w=sqrt((ofs-worst).^2);
% w=abs(ofs-worst);
w=w./(max(w)+eps);
