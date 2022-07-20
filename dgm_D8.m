function [y, x] = dgm_D8(sqW, sq_grips)
%% D8 routing 
%read in: sqW 3x3 elev grid and sq_grips 3x3 binary drain grid 1=drain 0=no drain 
%read out: x and y which are the index values for steepest downslope cell in 3x3 grid

W_grads= -1*(sqW-sqW(2,2)); %elev differences for each cell in sq +ve downslope

%the divide the diagonals by sqrt(2) so that they have correct gradient(pythag NB this is a scaling)
W_grads(1,1)=W_grads(1,1)/sqrt(2); 
W_grads(1,3)=W_grads(1,3)/sqrt(2);
W_grads(3,1)=W_grads(3,1)/sqrt(2);
W_grads(3,3)=W_grads(3,3)/sqrt(2);

%normalise slopes
scaler=max(max(W_grads)); %find max elev diff
sq_grads=W_grads/scaler; %scale so that max elev diff is 1
sq_weights=zeros(3,3); % 1st step for gradient weighting it is zeros
sq_weights(sq_grads>0)=(sq_grads(sq_grads>0)); %logical indexing for cells where sq_grads>0 replace the zeros is sq_weights with the sq_grads value
sq_weights2=sq_weights.*sq_grips; %NEW CODE sq_weights2 MAKES ALL THE NON GRIP CELLS ZERO
if sum(sum(sq_weights2))>0
    sq_weights=sq_weights2; 
end %.*sq_grips; %MAKES ALL THE NON GRIP CELLS ZERO
[y, x]=find(sq_weights==max(max(sq_weights)),1,'first'); %find cell ref of steepest downslope cell