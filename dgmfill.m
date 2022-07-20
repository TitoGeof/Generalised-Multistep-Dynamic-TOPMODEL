
function [W]=dgmfill(Z)

% %%New input
% clear all
% %if exist('gdem','var') ~= 1
% gdem=dlmread('oifdtm.asc',' ',7,1); % elevation [m]
% %i=find(gdem==-99999); gdem(i)=NaN;
% %end
% 
% Z = gdem(300:600,200:800);

%% Nick's code

Z_size=size(Z);
rows=Z_size(1,1);
cols=Z_size(1,2);

W=Z;

use_fill=0;% '1' means *RESET* the inner cells of W to the filler amount, otherwise
% the cells are reset to their source values + the pit tolerance - see below.

maxit=max(max(Z))+1;

W(2:rows-1,2:cols-1)=W(2:rows-1,2:cols-1)+maxit;

Zcol=reshape(Z,rows*cols,1);
Wcol=reshape(W,rows*cols,1);
diffs=Wcol-Zcol;
old_diffs=std(diffs);
Wstart=W;
Wstart_col=Wcol;
Wlast_col=Wcol;
clear Zcol Wcol diffs 

goes=100;%20%20%00
over_goes=33;
epsilon=0.001;

change_array=zeros(rows,cols);
trig=1;
fip_list=[1 6 3 8 2 7 4 5]';
%fip_list=[9 14 12 15 16 11 13 10]';

%fip_list=[13 14  2 7 4 5 15 12 1 9 11 6 3 16 8 10]';
fipple_list=zeros(8,7);
ortho1=[1 2 cols-1 1 2 rows-1 1];
ortho2=[2 2 cols-1 1 rows-1 2 -1];
ortho3=[3 cols-1 2 -1 2 rows-1 1];
ortho4=[4 cols-1 2 -1 rows-1 2 -1];
ortho5=[5 2 rows-1 1 2 cols-1 1];
ortho6=[6 rows-1 2 -1 2 cols-1 1];
ortho7=[7 2 rows-1 1 cols-1 2 -1];
ortho8=[8 rows-1 2 -1 cols-1 2 -1];

fipple_list(1,:)=ortho1;
fipple_list(2,:)=ortho6;
fipple_list(3,:)=ortho3;
fipple_list(4,:)=ortho8;
fipple_list(5,:)=ortho2;
fipple_list(6,:)=ortho7;
fipple_list(7,:)=ortho4;
fipple_list(8,:)=ortho5;

go_check=0;

for go=1:1:goes
    changed=0;
    counter=mod(go-1,8)+1;
    orth_num=fipple_list(counter,1);
    a1=fipple_list(counter,2);
    a2=fipple_list(counter,3);
    a_gap=fipple_list(counter,4);
    b1=fipple_list(counter,5);
    b2=fipple_list(counter,6);
    b_gap=fipple_list(counter,7);
    
    if orth_num<5
       swap_ij=1;
    else swap_ij=0;
    end
        
    
    for i=a1:a_gap:a2
        for j=b1:b_gap:b2
            if swap_ij==0
                ii=i;
                jj=j;
            else  ii=j;
                jj=i;
            end
            Zc=Z(ii,jj);
            Wc=W(ii,jj);
            sq=W(ii-1:ii+1,jj-1:jj+1);
            sq2=sq;
            sq2(2,2)=maxit;
            if Zc>=min(min(sq2))+epsilon
                W(ii,jj)=Zc;
                changed=changed+1;
                %change_array(i,j)=1;
            else if Wc>min(min(sq2))+epsilon
                    W(ii,jj)=min(min(sq2))+epsilon;
                end
            end
        end       
    end
    Wcol=reshape(W,rows*cols,1);
    diffs=Wlast_col-Wcol;
    std_diffs=std(diffs);
    
    if std_diffs==old_diffs
       ([ 'Breaking goes at go ' int2str(go) ])%', overgo ' int2str(over_go) '.' ])
       break
    else Wlast_col=Wcol;
         old_diffs=std_diffs;
    end
    %fred=([ 'W' int2str(over_go) ]);
    %eval([ fred '=W;' ]);
end