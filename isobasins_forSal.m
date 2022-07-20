function [out,DB]=isobasins_forSal(FD,Athresh,censor_frac)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Isobasins %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This is a Matlab implementation of the Whitebox Isobasins function but is
%'inspired by' rather than a direct translation. It generates a set of
%isobasins (catchments of user-defined area), some are first order
%catchments (with no upslope contributors) others have multiple upslope
%contributing catchments. This code uses functions from Topotoolbox v2. The
%current code is slow to run (ca 10 mins for 10^6 cells) but does appear
%to reproduce isobasins fairly consistently.

%% Inputs:
%FD, the flow direction object from topotoolbox v2
%Athresh, the desired average area for the isobasins (in cells at the mo)
%censor_frac, the fraction of the upslope area to censor so that hillslope cells aren't examined as target cells (suggest 0.5) lower=faster but can result in incomplete basin identification
%% Outputs:
%out, the outlet locations as a grid with dimensions of FD, cells have the id of the basin they are the outlet for or zero if not an outlet

%% Calculate upslope area and censor lowest areas.
A  = flowacc(FD);
A(A<(Athresh.*censor_frac))=0;

%% Initialise
out=zeros(size(A));
outletID=1;
nind=length(FD.ix);
for i=1:nind
    ind=FD.ix(i); % get index of target cell. NOTE:ix is the topologically sorted node list of givers, ixc is the same for receivers (in this single flowpath setup)
    if not(isnan(A(ind))) %don't evaluate noData cells
        if A(ind)>0 %don't evaluate cells with area set to zero in the censoring step
            flag=false; %initialise the flag that we use to jump out of while loop
            my_ind=ind; %initialise index that will march down the flowpath
            while not(flag)
                %% March down the flowpath from the target cell
                if sum(FD.ix==my_ind)==1
                    my_ind=FD.ixc(FD.ix==my_ind); % if there is a downslope neighbour then move to it. Note: ix is the giver, ixc is the receiver
                else
                    flag=true; %if no upslope neighbours jump out of while loop
                end
               
                %% Stop when the isobasin theshold is exceeded
                if A(my_ind)>=Athresh   %
                    %if the catchment area (UCA) of the cell greater than the desired isobasin area then get the id and UCA of the upslope neighbour with the largest UCA
                    icList=FD.ix(FD.ixc==my_ind); %find all the upslope cells
                    AsList=A(icList); %find their UCA
                    MaxAs=max(AsList); %find the max of the UCA
                    Maxic=icList(AsList==MaxAs); %get the id of the max uca cell
                    Maxic=Maxic(1); %take the first of these values
                   
                    %% Jump out of the while loop if any upslope neighbour cell has UCA above threshold
                    %because if this is the case you aren't on the dominant
                    %flowpath (i.e. in the river...) so you should try again
                    if MaxAs>Athresh
                        flag=true;
                    else
                        %% Choose the cell with area closest to isobasin threshold                        
                        d1=abs(A(my_ind)-Athresh);    %this is the cell you stopped marching down at
                        d2=abs(MaxAs-Athresh);          %this is the largest upslope neighbour
                        if d1>d2, out_ind=Maxic;
                        else, out_ind=my_ind;
                        end
                        out(out_ind)=outletID; %write label the outlet point in the grid
                       
                        %reduce the UCA along its entire flowpath to the domain edge by the isobasin area.
                        [myIX,~] = flowpathextract(FD,out_ind);     %find indices of downstream cells
                        A(myIX)=A(myIX)-Athresh;                %reduce their UCA by the isobasin area
                       
                        outletID=outletID+1; %increment outlet id
                        A(A<0)=0; %avoid negative areas
                    end
                end %if the area theshold is exceeded
            end %while to stop march down flowpath            
        end %if censoring cells with small UCA
    end %if censoring no-data cells
   
%     if mod(0.1*nind,i)
%         done=i/nind*100 %a progress message
%     end
   
end %for loop through the domain in topological order using ix in FD.


%% give the outlets map to drainage basins in topotoolbox and that will assign each cell in the grid to its associated basin (i.e. outlet)
IX=find(out(:)>0);              %get a list of indices for the isobasin outlets
DB = drainagebasins(FD,IX);     %calculate drainage basins for these outlets

end %end of function