function [allFiltered,newCenters,mu,sendPre,sendPost,discardedN,discardedRPre,discardedRPost,ribbon]=analyzeData(nuclei,presynaptic,postsynaptic,threshold,medRange,startValue,stopValue,sensitivity,rangeN,rangeR,radius,epsilon,minGroup,voxel,dimensions,allData)

%%%%THIS IS JUST THE FUNCTION FORM OF LOADDATA. IT IS USED IN THE GUI FOR
%%%%DOING EVERYTHING IN LOAD DATA.

% dataFile='IB09MidControlLeft2.czi';
% 
% 
% nuclei=1;
% presynaptic=2;
% postsynaptic=3;

%  threshold=[.1,.08,.08,.05];
% medRange=[8,8;3,3;3,3;4,4];

% plainData=allData;

%filter data by channel

allFiltered=zeros(dimensions(1),dimensions(2),dimensions(3),'logical');
% medRange=[8,8;7,7;7,7;8,8];
for b =1:4
    [allFiltered(:,:,:,b)]=inputAndThreshold(threshold(b),medRange(b,:),allData(:,:,:,b),b==nuclei);
end
% startValue=1;
% stopValue=20;


% sensitivity=.972;
% range=[1200,1500];
% radius=[50,90];
%finds circles on each slice for hte nucleus channel using the filtered
%data
 [storeCenters,storeRadii]=viewPreliminaryData(allFiltered(:,:,:,nuclei),rangeN,sensitivity,stopValue,startValue,radius);

 %clusters them together and outputs the final nucleus centers and radii
  [newCenters,mu,discardedN]=clusterNuclei(storeCenters,storeRadii,voxel);

  

%      rangeR=[1000,2048];
    ribbonSlices=[presynaptic,postsynaptic];
%     epsilon=10;
%     minGroup=5;
    ribbon=struct([]);
    noFit=struct([]);
    
    
    %go through the presynaptic ribbon and postsynaptic density channels
    %for all slices using the filtered data
    for i = 1:2
        %Find all ribbon locations on each slice
        [ribbon(i).points]=ribbonStuff(allFiltered(:,:,:,ribbonSlices(i)),epsilon,minGroup,rangeR,startValue,stopValue);

        %Group the ribbon locations together to create a 3D ribbon location
        %array. This part does nto work very well yet.
        [ribbon(i).grouped,noFit]=ribbonAnalysis(ribbon(i).points);
        if i==1
            discardedRPre=noFit;
        else
            discardedRPost=noFit;
        end
    end


%Ignore this part. It is just used to save time when I am experimenting
%with the GUI
    
    
[sendPre,sendPost]=ClusteredToSlice(ribbon,dimensions);
    
    
    
    
%     toGraph=struct([]);
%     for i = 1:size(ribbon.a,1)
%     end
