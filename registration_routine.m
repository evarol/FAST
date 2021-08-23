clear all
clc
close all
addpath ./tiff_loading/jsonlab/
addpath ./tiff_loading/utilities/
addpath(genpath('./tiff_loading/Fiji.app'));
javaaddpath('./tiff_loading/Fiji.app/mij.jar')

minmax = @(x)((x-min(x(:)))./max(x(:)-min(x(:))));
minmaxnnz = @(x)(max((x-min(x(x>0))+1)./max(x(:)-min(x(x>0))+1),0));
% load register_this.mat
z=6;
sigma=2;

%% load data
if ~exist('starmapfile')
    disp('Select confocal image...');
    [starmapfile,starmappath] = uigetfile('*','Select file to load','MultiSelect','on');
    starmap=double(load_tiff([starmappath starmapfile]));
    starmap=double(starmap(:,:,:,2));
end

if ~exist('starmapannotationfile')
    disp('Select confocal image annotation (if available)...');
    [starmapannotationfile,starmapannotationpath] = uigetfile('*','Select file to load','MultiSelect','on');
    if ~isequal(starmapannotationfile,0);s_annotation=load([starmapannotationpath starmapannotationfile]);end
end

if ~exist('zstackfile')
    disp('Select zstack image...');
    [zstackfile,zstackpath] = uigetfile('*','Select file to load','MultiSelect','on');
    zstack=double(load_tiff([zstackpath zstackfile]));
end

if ~exist('zstackannotationfile')
    disp('Select zstack image annotation (if available)...');
    [zstackannotationfile,zstackannotationpath] = uigetfile('*','Select file to load','MultiSelect','on');
    if ~isequal(zstackannotationfile,0);z_annotation=load([zstackannotationpath zstackannotationfile]);end
end



starmap_scale = [1272.7922/2048 1272.7922/2048 44/22];
% zstack_scale = [1000/700 1000/700 492/164];
zstack_scale = [1000/700 1000/700 2];
starmap_size = size(starmap).*starmap_scale;
zstack_size = size(zstack).*zstack_scale;

starmap=max(imresize3(starmap,starmap_size),0);
zstack=max(imresize3(zstack,zstack_size),0);

if ~isequal(starmapannotationfile,0)
Spoints=s_annotation.points.*starmap_scale;
starmap_filt=starmap;
end

if ~isequal(zstackannotationfile,0)
Zpoints=z_annotation.points.*zstack_scale;
Zstack_filt=zstack;
end


if isequal(starmapannotationfile,0)
%% pretty happy with this starmap segmentation
Istarmap=log1p(starmap);
Bstarmap=zeros(size(Istarmap));
for t=1:size(Istarmap,3)
    It=imtophat(Istarmap(:,:,t),strel('disk',11));
    starmap_filt(:,:,t)=It;
    Itg=imgaussfilt(It,sigma);
    Itgz=double(zscore(Itg,[],'all')>z);
    figure(99)
    imagesc([minmax(Istarmap(:,:,t)) minmax(It);minmax(Itg) minmax(Itgz)+minmax(It)]);
    Bstarmap(:,:,t)=Itgz;
    drawnow
end

Lstarmap=bwlabeln(Bstarmap);
starmapstats=regionprops3(Lstarmap);

Lfstarmap=Lstarmap;
[a,b]=ismember(Lstarmap(:),find(starmapstats.Volume<1200));
Lfstarmap(a==1)=0;
figure(1)
imagesc([minmax(max(log1p(starmap),[],3)) minmax(max(log1p(starmap),[],3)) + max(Lfstarmap>0,[],3)]);
Sstats = regionprops3(Lfstarmap);
Spoints = Sstats.Centroid;Svolumes=Sstats.Volume;Svolumes(any(isnan(Spoints),2),:)=[];Spoints(any(isnan(Spoints),2),:)=[];
end

if isequal(zstackannotationfile,0)
%% pretty okay with this segmentation too
Izstack=log1p(zstack);
Bzstack=zeros(size(Izstack));
for t=1:size(Izstack,3)
    It=log1p(exp(imtophat(log1p(zstack(:,:,t)),strel('disk',11)))-exp(imtophat(log1p(zstack(:,:,t)),strel('disk',4))));
    Zstack_filt(:,:,t)=It;
    Itgz=double(zscore(It,[],'all')>z);
    L=bwlabeln(Itgz);
    stats=regionprops3(L);
    [a,b]=ismember(L(:),find(stats.Volume<100));
    L(a==1)=0;
    Itgz=double(L>0);
    figure(99)
    imagesc([minmax(Izstack(:,:,t)) minmax(It);minmax(Itgz) minmax(Itgz)+minmax(It)]);title(num2str(t));
    Bzstack(:,:,t)=Itgz;
    drawnow
end

Lzstack=bwlabeln(Bzstack);
zstackstats=regionprops3(Lzstack);

Lfzstack=Lzstack;
[a,b]=ismember(Lzstack(:),find(zstackstats.Volume<1200));
Lfzstack(a==1)=0;
figure(2)
imagesc([minmax(max(Zstack_filt,[],3)) minmax(max(Zstack_filt,[],3)) + 0.5*max(Lfzstack>0,[],3)]);
Zstats = regionprops3(Lfzstack);
Zpoints = Zstats.Centroid;Zvolumes=Zstats.Volume;Zvolumes(any(isnan(Zpoints),2),:)=[];Zpoints(any(isnan(Zpoints),2),:)=[];
end
%% Point set registration
figure(3)
subplot(1,2,1)
hold on
plot3(Spoints(:,1),Spoints(:,2),Spoints(:,3),'b.','MarkerSize',20);axis equal;axis tight
plot3(Zpoints(:,1),Zpoints(:,2),Zpoints(:,3),'r.','MarkerSize',20);axis equal;axis tight
legend('Confocal','Zstack');
set(gca,'FontWeight','bold','FontSize',20,'TickLength',[0 0]);set(gcf,'Color','w');

[bhat,P,inlier_set,q]=rrwoc_wahba(Spoints,Zpoints,0.9,30,1,2000);

% [vol_dist,P]=volumetric_distance_2(Spoints,Zpoints,10);

% [vol_dist]=knn_distance(Spoints,Zpoints,50,30);
% M=munkres(-vol_dist);
% 
% [R,T,beta]=weighted_wahba(Spoints,M*Zpoints,max(M.*exp(vol_dist./2*0.1^2),[],2));
% 
% bhat=[R;T];

Zpointshat=[Spoints ones(size(Spoints,1),1)]*bhat;

[P.matched_cells_y,P.matched_cells_x]=find(pdist2(Zpoints,Zpointshat)<= 30);
figure(3)
subplot(1,2,2)
hold on
plot3(Zpointshat(:,1),Zpointshat(:,2),Zpointshat(:,3),'b.','MarkerSize',20);axis equal;axis tight
plot3(Zpoints(:,1),Zpoints(:,2),Zpoints(:,3),'r.','MarkerSize',20);axis equal;axis tight
legend('Registered Confocal','Zstack');
set(gca,'FontWeight','bold','FontSize',20,'TickLength',[0 0]);set(gcf,'Color','w');


starmap_histmatched=linhistmatch(log1p(starmap_filt),minmax(max(log1p(Zstack_filt),[],3)),200,'regular');
tic
tform = affine3d([[bhat(1:3,1:3);bhat(4,:)] [0;0;0;1]]);
starmap_moved=imwarp(starmap_histmatched,tform,'Outputview',imref3d(size(Zstack_filt)));
toc


figure(4)

A=max(starmap_moved,[],3);
B=max(minmax(log1p(Zstack_filt)),[],3);

I(:,:,1)=[A A.*B zeros(size(B))];
I(:,:,2)=[zeros(size(A)) zeros(size(A)) zeros(size(B))];
I(:,:,3)=[zeros(size(A)) A.*B B];

imagesc(I)
hold on
plot(Zpointshat(Zpointshat(:,1)<size(A,1),1),Zpointshat(Zpointshat(:,1)<size(A,1),2),'y.','MarkerSize',10)
plot(Zpointshat(P.matched_cells_x,1)+size(A,2),Zpointshat(P.matched_cells_x,2),'y.','MarkerSize',10);
for i=1:length(P.matched_cells_y)
 plot([Zpointshat(P.matched_cells_x(i),1)+size(A,2) Zpoints(P.matched_cells_y(i),1)+size(A,2)],[Zpointshat(P.matched_cells_x(i),2) Zpoints(P.matched_cells_y(i),2)],'w-','LineWidth',2);   
end
plot(Zpoints(P.matched_cells_y,1)+size(A,2),Zpoints(P.matched_cells_y,2),'c.','MarkerSize',10);
plot(Zpoints(:,1)+size(A,2)+size(A,2),Zpoints(:,2),'c.','MarkerSize',10);