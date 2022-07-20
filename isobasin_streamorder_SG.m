function [B,NB]=isobasin_streamorder_SG(DEM,FDs,CHthresh)
% [R,ch,~]=channel_reaches(DEM,FDs,CHthresh);
% %define a channel network
% ch       = 0*DEM.Z;
% ch(As*cs^2>=CHthresh)=1;
% [~,ch]=channel_reaches(DEM,FDs,CHthresh);
% %obtain stream order
% S               = streamorder(FDs,ch);
% %retain unique IDs
% B0              = drainagebasins(FDs,S,1) ;
% B               = B0+1;
% %single flow direction for surface[object]
% FDs                 = FLOWobj(DEM);
%get stream object
S               = STREAMobj(FDs,'minarea',CHthresh);
%split at intersections
S               = split(S);
B0              = drainagebasins(FDs,S) ;
B               = B0+1;
B(isnan(DEM.Z)) = NaN;
%adjust basin ID to go from 1 to Nr (Nr would be number of channel reaches too)
NB          = unique(B);
NB(isnan(NB) | NB==0)=[];
for ii=1:length(NB)
    cond=B == NB(ii);
    B(cond) = ii;
    NB(ii) = ii;
end
%convert class to int32
Y  = typecast(B(:),'int32');
NB = typecast(NB,'int32');
%reshape into matrix
Y  = reshape(Y,(size(B)));
%convert into double
B  = double(Y);
NB = double(NB);
NB = max(NB);