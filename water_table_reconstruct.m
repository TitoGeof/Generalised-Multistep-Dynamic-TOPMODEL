function DEM=water_table_reconstruct(DEM,As,FD,alpha,datum)
%--------------------------------------------------------------------------
%reconstructs water table from subsurface hydraulic gradient 'alpha [degree]'
%--------------------------------------------------------------------------
% DEM   : surface DEM object [GridOBJ]
% FD    : single flow direction object [FlowOBJ]
% As    : single direction flow accumulation
% alpha : subsurface hydraulic gradient (from Hjredt eta al. [2004])
% datum : elevation difference between surface and water table @outlet-cell 
%--------------------------------------------------------------------------
%giver cell IDs (flip to start from outletcell)
giv         = flipud(FD.ix);
%receiver cell IDs(flip to start from outletcell)
rec         = flipud(FD.ixc);
%distance between givers and their receivers (either cs or sqrt(2)*cs)
dx          = getdistance(giv,rec,DEM.size,DEM.cellsize);
%generate a vector of unique receivers to loop through (keep order intact)
recu        = unique(rec,'stable');
% recu=rec;
%number of cells to loop through
nCells      = length(recu);
%initialse subsurface DEM (sparse)
% Z=spalloc(DEM.size(1),DEM.size(2),not(isnan(nnz(DEM.Z))));
Z=DEM.Z*0;
%find outlet cell
outletID    = find(As==max(As(:)));
%set outlet elevation equal to surface - datum
Z(outletID) = DEM.Z(outletID) - datum;
%loop 
for ii      = 1:nCells
    %find the givers to this receiver
    ID         = find(rec==recu(ii)); 
    %back calculate giver cell elevations
    Z(giv(ID)) = Z(recu(ii)) + dx(ID).*tand(alpha(recu(ii)));
end
Z=full(Z)+DEM.Z*0;
class(Z)
%subsurface DEM gridOBJ
DEM.Z     = Z;


% %giver cell IDs (flip to start from outletcell)
% giv         = flipud(FD.ix);
% %receiver cell IDs(flip to start from outletcell)
% rec         = flipud(FD.ixc);
% %distance between givers and their receivers (either cs or sqrt(2)*cs)
% dx          = getdistance(giv,rec,DEM.size,DEM.cellsize);
% %generate a vector of unique receivers to loop through (keep order intact)
% recu        = unique(rec,'stable');
% % recu=rec;
% %number of cells to loop through
% nCells      = length(recu);
% %initialse subsurface DEM (sparse)
% Z=spalloc(DEM.size(1),DEM.size(2),not(isnan(nnz(DEM.Z))));
% %find outlet cell
% outletID    = find(As==max(As(:)));
% %set outlet elevation equal to surface - datum
% Z(outletID) = DEM.Z(outletID) - datum;
% %loop 
% for ii      = 1:nCells
%     %find the givers to this receiver
%     ID         = rec==recu(ii);
%     %back calculate giver cell elevations
%     Z(giv(ID)) = Z(recu(ii)) + dx(ID).*tand(alpha(recu(ii)));
% end
% Z=full(Z)+DEM.Z*0;
% % DEM.Z(not(isnan(nnz(DEM.Z))))=Z;
% %subsurface DEM gridOBJ
% DEM.Z     = Z;