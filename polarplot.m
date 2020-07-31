function [] = polarplot(file_nuc,file_kev,file_ido)
% outputs polarplot comparing files: file_kev and file_ido and defining # cluster
% with file_nuc
% requires k_med or kcluster function to be added to path for polarplot to
% work
%% rotate segmented single nucleus and its respective ribbons to be perpindicular to xy plane
% rotation of unique axis
%% decide which thresholding/clustering method to use
num = 'Which clustering method xaxis or kcluster or kmed ? '; 
method = input(num,'s');
%% load data  nucleus kevins
kxyz=readtable(file_nuc) ;
Kxyz=table2array(kxyz);
KXYZ=Kxyz(:,1:3) ;   
nucx=(KXYZ(:,1));
nucy=(KXYZ(:,2));
nucz=(KXYZ(:,3));
nucx=sort(nucx);
nucy=sort(nucy);
nucz=sort(nucz);
%% load data presynaptic ribbons
kxyz =readtable(file_kev) ;
ixyz=readtable(file_ido) ;
 
Kxyz=table2array(kxyz);
KXYZ=Kxyz(:,1:3) ;   
kx=(KXYZ(:,1));
ky=(KXYZ(:,2));
kz=(KXYZ(:,3));

Ixyz=table2array(ixyz(:,1:3));
IXYZ=Ixyz(:,1:3) ;   
ix=(IXYZ(:,1));
iy=(IXYZ(:,2));
iz=(IXYZ(:,3));
%% segment ribbons to their nucleus
%create struct to partition presynaptic ribbons to groups belonging to a nucleus
for t=1:2 % repeat code for ido and kevin data
   
      if t==2 % to switch to  idos data and reduce repetive code lines
        keven=[kx ky kz]; % save kevins data
        kx=ix;
        ky=iy;
        kz=iz;
      end
KX=struct([]);
KY=struct([]);
KZ=struct([]);
thresh=[];
a=sort(nucx); % rearrange from least to greatest top to bottom
for i=1:length(nucx)-1
thresh(i)=(a(i+1)+a(i))/2;
end
thresh=[0 thresh];

for k=1:length(thresh)-1
    kx1=[];
    ky1=[];
    kz1=[];
    
    kx2=[];
    ky2=[];
    kz2=[];
for j=1:length(kx)
    if kx(j)<thresh(k+1)&& kx(j)>thresh(k)
        tempx=kx(j);
        tempy=ky(j);
        tempz=kz(j);
        kx1=[kx1 tempx];
        ky1=[ky1 tempy];
        kz1=[kz1 tempz];
    elseif kx(j)>max(thresh)
        tempx=kx(j);
        tempy=ky(j);
        tempz=kz(j);
        kx2=[kx2 tempx];
        ky2=[ky2 tempy];
        kz2=[kz2 tempz];
    end
    end
    KX{k}=kx1;
    KY{k}=ky1;
    KZ{k}=kz1;
end
 KX{k+1}=kx2;
 KY{k+1}=ky2;
 KZ{k+1}=kz2;
 %% kmdedoids clustering to replace x axis threshhold
 switch method
     case 'kmed'
        [KX KY KZ Centroid]=kmed(kx,ky,kz,length(nucx)); %same issue
 %% k clustering to replace x -axis threshhold
 % issue withk clusteringn is that idos and kevins cluster differently
 %so cannot compare since the clustering is not similar to both
     case 'kcluster'
        [KX KY KZ Centroid]=kcluster(kx,ky,kz,length(nucx));
 end
 %% find centroid of all clusters
 for i=1:length(KX)
 avgKX(i)=mean(KX{i});
 avgKY(i)=mean(KY{i});
 avgKZ(i)=mean(KZ{i});
 end
% avgKX=Centroid(:,1);
 %avgKY=Centroid(:,2);
 %avgKZ=Centroid(:,3);
%% plot segmented presynaptic nucleus to its respecitve nucleus
if t==1 % to avoid repetiveness while I run ido data
num = 'how many nucleus do you want plotted? ';
number = input(num);
end
figure (1)
for i=1:number
       switch t
           case 1
               prompt = 'What nucleus do you want to isolate in Kevin? ';
           case 2
               prompt = 'What nucleus do you want to isolate in Ido? ';
       end
nuc = input(prompt);
switch t
    case 1
        scatter3(KX{nuc},KY{nuc},KZ{nuc},'ro')
        hold on
        scatter3(nucx(nuc),nucy(nuc),nucz(nuc),'ko')
        scatter3(avgKX(nuc),avgKY(nuc),avgKZ(nuc),'k*') % centroid of specifc ribbon cluster
    case 2
        hold on
          scatter3(nucx(nuc),nucy(nuc),nucz(nuc),'ko')
          scatter3(avgKX(nuc),avgKY(nuc),avgKZ(nuc),'k*')
          if i==max(number)
              legend('Kevin','Nucleus','Centroid')
              legend('boxoff')
          end
          scatter3(KX{nuc},KY{nuc},KZ{nuc},'DisplayName','Ido','MarkerEdgeColor','b') % ido data
end
end
hold off
%xlim([-100,100])
%ylim([-100,100])
zlim([-20,40])
xlabel('x')
ylabel('y')
zlabel('z')
%% transition centroid and its nucleus to center (0,0)
if number==1
    figure (2)
    transKX=avgKX(nuc)-nucx(nuc);
      transKY=avgKY(nuc)-nucy(nuc);
      transKZ=avgKZ(nuc)-nucz(nuc);
      switch t
          case 1
        scatter3(transKX,transKY,transKZ,'r*')
        hold on
        scatter3(0,0,0,'ko')
        xlabel('x')
        ylabel('y')
        zlabel('z')
        xlim([-20,20])
        ylim([-20,20])
        zlim([-20,20])
          case 2
              hold on
              scatter3(transKX,transKY,transKZ,'b*')
              legend('kevin','nucleus','Ido')
      end
      hold off
end
 %% First dot product angle between two lines to find two angles for rotation
% goal is to average the single cluster to find centroid of cluster and
% make the vector of this centroid-nucleus line parallel with z-axis perp. to xzplane
if number==1
a=[transKX 0 transKZ ];% centroid XZ vector projection
b=[0 0 1];% z-axis
XZangle=acosd(dot(a,b)/(norm(a)*norm(b)));

%error code below since XZ rotation makes YZ angle ~= starting position after
%first rotation so need to do dot product after first rotation to figure out
%correct angle.
%{ 
a=[0 transKY transKZ ];% centroid YZ vector projections
YZangle=acosd(dot(a,b)/(norm(a)*norm(b)));
%}
%% rotate centroid about y-axis

%Around Y-axis: 
theta=(XZangle); % degree only
X = transKX*cosd(theta) + transKZ*sind(theta);
Y = transKY;
Z = transKZ*cosd(theta) -transKX*sind(theta);
%% second angle dot product
a=[0 Y Z ];% centroid YZ vector projections
YZangle=acosd(dot(a,b)/(norm(a)*norm(b)));
%% rotate centroid about x -axis
%Around X-axis: 
temp=[X; Y ;Z];
theta=(YZangle);% degree
X = temp(1,:);
Y =temp(2,:)*cosd(theta) - temp(3,:)*sind(theta);
Z = temp(2,:)*sind(theta) + temp(3,:)*cosd(theta);
end
%% plot to verify that centroid vector is parallel to z axis
if number==1
    figure (3)
    switch t
        case 1
            scatter3(X,Y,Z,'r*')
            cx=X+nucx(nuc); % save new centroid coord. and transition back
            cy=Y+nucy(nuc);
            cz=Z+nucx(nuc);
            hold on
        
            scatter3(0,0,0,'ko')
            xlabel('x')
            ylabel('y')
            zlabel('z')
            xlim([-20,20])
            ylim([-20,20])
            zlim([-20,20])
        case 2
            hold on
            scatter3(X,Y,Z,'b*')
            legend('kevin','nucleus','Ido')
    end  
    hold off
end
%% plot desired segmented correctly oriented cluster and nucleus polarplot
%% transition desired cluster to center around desired nucleus
if number==1
    transKX=KX{nuc}-nucx(nuc);
      transKY=KY{nuc}-nucy(nuc);
      transKZ=KZ{nuc}-nucz(nuc);
%%
%Around Y-axis: 
theta=(XZangle); % degree only

X = transKX*cosd(theta) + transKZ*sind(theta);
Y = transKY;
Z = transKZ*cosd(theta) -transKX*sind(theta);

temp=[X; Y ;Z];
%Around X-axis:
theta=(YZangle); % degree only
X = temp(1,:);
Y =temp(2,:)*cosd(theta) - temp(3,:)*sind(theta);
Z = temp(2,:)*sind(theta) + temp(3,:)*cosd(theta);
%% transition cluster back to original position
 X=X+nucx(nuc);
 Y=Y+nucy(nuc);
 Z=Z+nucz(nuc);
%% final plot verification and final polarplot of nucleus and ribbon of interest
figure(4)
switch t
    case 1
        scatter3(X,Y,Z,'ro')
        hold on
        scatter3(nucx(nuc),nucy(nuc),nucz(nuc),'ko')
        scatter3(avgKX(nuc),avgKY(nuc),avgKZ(nuc),'k*') % centroid of specifc ribbon cluster
        scatter3(cx,cy,cz,'k*') % new location of centroid of specifc ribbon cluster
        zlim([-20,40])
        xlabel('x(um)')
        ylabel('y(um)')
        zlabel('z(um)')
    case 2
        hold on
        legend('kevin','Nucleus','centroid')
        scatter3(X,Y,Z,'DisplayName','Ido','MarkerEdgeColor','b')
end
hold off
r=(sqrt((X-nucx(nuc)).^2+(Y-nucy(nuc)).^2))'; % center ribbons to nucleus
phi=(atand((Y-nucy(nuc))./(X-nucx(nuc))).*(2*pi/360))';% converts to radians

for i=1:length(KX{nuc})
    if X(i)-nucx(nuc)<0
       phi(i)=phi(i)+pi;% fixes the issue of atand limited to range -90 to 90
    end
end
figure(5)
switch t
    case 1
        %polarplot(phi,r,'ro')% accepts phi as radians and converts to degrees
        mypolar(phi,r,'ro')
    case 2
        hold on
        %polarplot(phi,r,'bo')% accepts phi as radians and converts to degrees
        mypolar(phi,r,'bo')
end
hold off
title (sprintf('Dataset %d Ribbons centered around Nucleus (%1fum,%1fum)',dataset,round(nucx(nuc),1),round(nucy(nuc),1)))
legend('Kevin','Ido')
legend('boxoff')
end
end
if t==2
kx=keven(:,1); % restore kevins values
ky=keven(:,2);% restore kevins values

kz=keven(:,3);% restore kevins values
end

end

