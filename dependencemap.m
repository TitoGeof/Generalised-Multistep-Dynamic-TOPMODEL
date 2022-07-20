function OUT = dependencemap(FD,varargin)





% if nargin == 2;
% SEED pixels are either supplied as logical matrix, GRIDobj, or linear
% index
SEED = varargin{1};
% isGRIDobj = isa(SEED,'GRIDobj');

% SEED is supposed to be supplied as linear index in the GRIDobj
ix   = varargin{1};
SEED = false(FD.size);
ix   = round(ix);
% if any(ix <= 0 | ix > prod(FD.size))
%     error('TopoToolbox:WrongInput',...
%         ['Linear indices must not be less or equal to zero \n' ...
%         'or larger than ' double2str(prod(FD.size)) '.']);
% end

SEED(ix) = true;

% if ~(exist(['dependencemap_mex.' mexext],'file') == 3)
% m implementation
ixtemp  = FD.ix;
ixctemp = FD.ixc;
for r = numel(ixtemp):-1:1
    SEED(ixtemp(r)) = SEED(ixtemp(r)) || SEED(ixctemp(r));
end

OUT=SEED;
% %% Prepare Output
% % empty GRIDobj
% OUT = copy2GRIDobj(FD);
% % write output to GRIDobj
% OUT.Z = SEED;
% OUT.zunit = 'logical';
% OUT.name  = 'dependence map';


end