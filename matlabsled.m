% This code is Victor Granlunds second best entry (Turbopulkka 2000) in
% Reaktors Traveling Santa challenge. Written in MATLAB R2020b. Note that
% it doesn't calculate the final distance, plug in your own script for
% visualization and distance calculation. Distance matrix sections might
% take 20 minutes to run, but the algorithm shouldn't take more than a minute
% to run. Make sure you have nicelist.txt in the same dir.
Korvatunturi=[1 68.073611 29.315278 0]; %Korvatunturi is set as row 1, because it allows us to have distance to and from korvatunturi for each child in the distance matrix.
Nicelist=table2array(readtable('nicelist.txt'));    %the array is filled with the rest of the children
T=[Korvatunturi;table2array(readtable('nicelist.txt'))];
 %% Generate distance matrix (Using a distance matrix is required so that the rest of the code can be run fast)
 T_distmatrix=zeros(height(T),height(T));
 for from = 1:height(T)
     for to = 1:height(T)
         T_distmatrix(from,to)=acosd(sind(T(to,2))*sind(T(from,2))+cosd(T(to,2))*cosd(T(from,2))*cosd(-T(from,3)+T(to,3)));
     end
     from %Tracking progress
 end
 % Put NaN on diagonal (So that minimization works without using exceptions)
 for diagonal = 1:height(T) %this does take a suprising amount of time, for:loops are super inefficient.
     T_distmatrix(diagonal,diagonal)=NaN;
     diagonal
 end
 % The distance matrix is large and takes some time to generate. It's useful to store a backup incase this section is run again by mistake.
 Backupdistancematrix=T_distmatrix; %you can also save it as a file once it's generated.
 %% Use backup to get T_distmatrix
 T_distmatrix=Backupdistancematrix;
%% Algorithm: Once the distance matrix has been generated, make sure to run this section only
        %Optimization variables. The given values are far from optimal.
maxdistance=            20;      %Controls the max distance to search for additional deliveries when the sack is close to full. Note that there is interplay with homebias, which replaces distance with a biased distance to neighbors.
maxdistanceexponent=    2;     %Controls how fast maxdistance comes into play when sack weight nears sackcapacity. The intuition behind this strategy is that if you have room for lots of packages more, it's fine to go a long way for the next package.
homebiasexponent=       2;    %The bias towards korvatunturi is 1 when the sack is full. This exponent controls how fast or slow it goes to 1.
% note that distances are in degrees, not kilometers
        
    %Code starts
ItemsLeft=Nicelist(:,1);
SackCapacity=10000*1000; %grams
clear Delivered     
DeliverRow=1;       %Tracker used in while loop (same rows as in final output)
while length(ItemsLeft)>0
    SackItems=[];
    SackWeight=0;
    PositionIndex=1;    %on index 1 in distance matrix is korvatunturi
    FirstSackItem=1;    %The item per run is chosen by farthest-neighbor, not nearest-neighbor
    while (SackWeight <= SackCapacity) & not(isempty(ItemsLeft))
        Candidates=zeros(length(ItemsLeft),3);  %preallocation
        for checker = 1:length(ItemsLeft)       %This loop creates a list with deliverable items, their biased distances, and masses.
            Candidates(checker,:)=[ItemsLeft(checker) T_distmatrix(PositionIndex,ItemsLeft(checker))-(((SackWeight-16000)/(SackCapacity-16000)))^homebiasexponent*(T_distmatrix(PositionIndex,1)-T_distmatrix(ItemsLeft(checker),1)) T(ItemsLeft(checker),4)];
            % Index     Biased-Distance    Mass. Note that the bias is negative when sack is filled less than 16000 grams. This means thesled will prefer packages that take you farther away from korvatunturi.
        end
        Candidates(isnan(Candidates))=0; %handles a special case
        if FirstSackItem == 1            %farthest-neighbor
            Candidates=sortrows(Candidates,2,'descend');
            FirstSackItem=0;
        else                             %nearest-neighbor
            Candidates=sortrows(Candidates,2,'ascend');
        end
        if min(Candidates(:,3)) > (SackCapacity - SackWeight)   %if none of the items fit in the sack, break now
            break
        end
        didntaddone = 1;
        for index = 1:size(Candidates,1)    %goes through all candidates to find a packages to add
            if (Candidates(index,3) <= (SackCapacity - SackWeight))&(((SackWeight/SackCapacity)^maxdistanceexponent*Candidates(index,2)) < maxdistance)
                SackItems=[SackItems Candidates(index,1)];
                SackWeight=SackWeight+Candidates(index,3);
                PositionIndex=Candidates(index,1);
                ItemsLeft = ItemsLeft(ItemsLeft~=Candidates(index,1));
                didntaddone=0;
                break           %if something has been added, we need to recalc distances.
            end
        end
        if didntaddone == 1
            break           %if no packages were acceptable to add to sack, stop trying to fill the sack and go to next run.
        end
    end
    if isempty(SackItems)   %if the sack filling loop finished without adding anything, stop trying to do deliveries. this avoids an unending loop in case of weird errors, and the results will be invalid.
        'invalid result, no candidates were deemed acceptable during sack filling.'
        break
    end
    Delivered{DeliverRow}=SackItems;
    DeliverRow=DeliverRow+1;
end
disp(DeliverRow-1) %display N of rows, i.e. how many runs from korvatunturi were made.
strdelivery=''; %Output formatting. Note that reaktors site actually wants txts, not csv:s. Just rename it before delivery or to what format your output processing script needs.
for a = 1:length(Delivered)
    if not (a==1)
        strdelivery=append(strdelivery,newline);
    end
    for b = 1:length(Delivered{a})
        if not (b == 1)
            strdelivery=append(strdelivery,'; ');
        else
            firstentry=0;
        end
        strdelivery=append(strdelivery,num2str(Delivered{a}(b)));
    end
end
writematrix(strdelivery,'matlabturbosledlist.csv');
%system('distance-calculator.py');

%This my second best code for the challenge, good enough to get on the leaderboard but
%not in the top ten. The code used for my entry victorg is a further
%improvement on this code, with about 20 lines added.