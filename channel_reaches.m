function [R,ch,NR]=channel_reaches(DEM,FDs,thresh)
%get stream object
S        = STREAMobj(FDs,'minarea',thresh);
%split at intersections
S        = split(S);
if not(isempty(S.ix))
  %obtain cell IDs of each reach
  reachIDs = S.orderednanlist;
  %Reach matrix
  R        = 0*DEM.Z;
  %grid  values
  I        = S.IXgrid;
  %seperate into cell arrays
  idx      = isnan(reachIDs);
  idy      = 1+cumsum(idx);
  idz      = 1:size(reachIDs,1);
  C        = accumarray(idy(~idx),idz(~idx),[],@(r){reachIDs(r,:)});
  %assign reach IDs in the matrix
  N        = size(C,1);
  for ii   = 1:N
    R(I(C{ii}))=N-ii+1;
  end
  %number of reaches
  Runq     = unique(R);
  NR       = nnz(not(isnan(Runq)));
  %overall channel map
  ch       = R>0;
else
  %if there's no channel with the specified threshold
  R        = zeros(size(DEM.Z));
  ch       = R;
  NR       = 0;
end