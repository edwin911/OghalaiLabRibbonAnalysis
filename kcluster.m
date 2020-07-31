function [KX,KY,KZ,Centroid] = kcluster(kx,ky,kz,nucleusCount)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% K clustering to improve threshholding
%% k cluster
%use kmeans clustering, setting the number of regions = to the number of
%found nuclei (number of rows in nuclei)
%iter=50;% number of iterations. higher the value the more accurate the k clustering
X = [kx ky kz];
KX=struct([]);
KY=struct([]);
KZ=struct([]);
iterations=200;
opts = statset('Display','final');
[idx,C] = kmeans(X,nucleusCount,'Distance','sqeuclidean',...
'Replicates',iterations,'Options',opts);
%'sqeuclidean' 'cityblock'
%{
for i=1:nucleusCount
Kx{i}=X(idx==i,1)';
Ky{i}=X(idx==i,2)';
Kz{i}=X(idx==i,3)';
end
KX =Kx;
KY=Ky;
KZ=Kz;
Centroid=C;
%}
%% reorganize data

tempx=struct([]);
tempy=struct([]);
tempz=struct([]);

kx=[];
for i=1:nucleusCount
    tempx{i}=X(idx==i,1);
    tempy{i}= X(idx==i,2);
    tempz{i}=X(idx==i,3);
    
    kx=[kx max(X(idx==i,1))];
end
KX=struct([]);
KY=struct([]);
KZ=struct([]);

for i=1:nucleusCount
    tempX=struct([]);
    tempY=struct([]);
    tempZ=struct([]);
    KX{i}=tempx{find(min(kx)==kx)}'
    KY{i}=tempy{find(min(kx)==kx)}'
    KZ{i}=tempz{find(min(kx)==kx)}'
    
    for j=1:length(tempx)
        if ~(max(tempx{find(min(kx)==kx)})==max(tempx{j}))
        tempX{j}=tempx{j};
        tempY{j}= tempy{j};
        tempZ{j}= tempz{j};
        else
            break
        end
    end
    for i=j:length(tempx)-1 
        tempX{i}=tempx{i+1};
        tempY{i}= tempy{i+1};
        tempZ{i}= tempz{i+1};
    end
    kx(find(min(kx)==kx))=[];
    tempx=tempX;
    tempy=tempY;
    tempz=tempZ;
end
Centroid=C;
end

