clear all
clc
close all

%% dependencies to load tiffs/nd2s
addpath ./tiff_loading/jsonlab/
addpath ./tiff_loading/utilities/
addpath(genpath('./tiff_loading/Fiji.app'));
javaaddpath('./tiff_loading/Fiji.app/mij.jar')


%% auxillary functions
minmax = @(x)((x-min(x(:)))./max(x(:)-min(x(:))));
minmaxnnz = @(x)(max((x-min(x(x>0))+1)./max(x(:)-min(x(x>0))+1),0));
vec=@(x)(x(:));

%% image parameters
confocal_scale = [1272.7922/2048 1272.7922/2048 44/22]; %x,y,z pixel dimension - user specified. VERY IMPORTANT to have real spatial dimensions
zstack_scale = [1000/700 1000/700 2]; %x,y,z pixel dimension - user specified. VERY IMPORTANT to have real spatial dimensions

%% Detection parameters
z=6; %z-score to distinguish bright foreground from background
sigma=2; %smoothing kernel to fill in gaps (to avoid spotty cells)

%% Registration parameters
radius_parameter = 200; % how big of a window around each cell we look for matching partners between modalities
tightness_parameter = 10; % how close do two registered cells need to be to be considered a match
max_depth_parameter = 50; % what is the maximum depth in z shall we look for matching cell neighbors
anglelimit=30; %maximum horizontal angle allowed

%% load data
if ~exist('confocalfile')
    disp('Select confocal image...');
    [confocalfile,confocalpath] = uigetfile('*','Select file to load','MultiSelect','on');
    confocal=double(load_tiff([confocalpath confocalfile]));
    confocal=double(confocal(:,:,:,2));
end

if ~exist('confocalannotationfile')
    disp('Select confocal image annotation (if available)...');
    [confocalannotationfile,confocalannotationpath] = uigetfile('*','Select file to load','MultiSelect','on');
    if ~isequal(confocalannotationfile,0);s_annotation=load([confocalannotationpath confocalannotationfile]);end
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


%% convert to micron dimensions
confocal_size = size(confocal).*confocal_scale;
zstack_size = size(zstack).*zstack_scale;

%% resize image such that 1 pixel = 1um^3
confocal=max(imresize3(confocal,confocal_size),0);
zstack=max(imresize3(zstack,zstack_size),0);

%% load annotations if available
if ~isequal(confocalannotationfile,0)
    Spoints=s_annotation.points.*confocal_scale;
    Iconfocal=log1p(confocal);
    for t=1:size(Iconfocal,3)
        It=imtophat(Iconfocal(:,:,t),strel('disk',11));
        confocal_filt(:,:,t)=It;
    end
end

%% load annotations if available
if ~isequal(zstackannotationfile,0)
    Zpoints=z_annotation.points.*zstack_scale;
    Izstack=log1p(zstack);
    for t=1:size(Izstack,3)
        It=log1p(exp(imtophat(log1p(zstack(:,:,t)),strel('disk',11)))-exp(imtophat(log1p(zstack(:,:,t)),strel('disk',4))));
        Zstack_filt(:,:,t)=It;
    end
end

%% Automatic cell detection
if isequal(confocalannotationfile,0)
    %% pretty happy with this confocal segmentation
    Iconfocal=log1p(confocal);
    Bconfocal=zeros(size(Iconfocal));
    for t=1:size(Iconfocal,3)
        It=imtophat(Iconfocal(:,:,t),strel('disk',11));
        confocal_filt(:,:,t)=It;
        Itg=imgaussfilt(It,sigma);
        Itgz=double(zscore(Itg,[],'all')>z);
        figure(99)
        imagesc([minmax(Iconfocal(:,:,t)) minmax(It);minmax(Itg) minmax(Itgz)+minmax(It)]);
        Bconfocal(:,:,t)=Itgz;
        drawnow
    end
    
    Lconfocal=bwlabeln(Bconfocal);
    confocalstats=regionprops3(Lconfocal);
    
    Lfconfocal=Lconfocal;
    [a,b]=ismember(Lconfocal(:),find(confocalstats.Volume<1200));
    Lfconfocal(a==1)=0;
    figure(1)
    imagesc([minmax(max(log1p(confocal),[],3)) minmax(max(log1p(confocal),[],3)) + max(Lfconfocal>0,[],3)]);
    Sstats = regionprops3(Lfconfocal);
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

%% filter out garbage for image based registration (not used right now)

% confocal_filt2=zeros(size(confocal_filt));
%
% for i=1:size(Spoints,1)
%     confocal_filt2(round(Spoints(i,2)),round(Spoints(i,1)),round(Spoints(i,3)))=100;
% end
% confocal_filt2=imgaussfilt3(confocal_filt2,5);
%
%
% Zstack_filt2=zeros(size(Zstack_filt));
%
% for i=1:size(Zpoints,1)
%     Zstack_filt2(round(Zpoints(i,2)),round(Zpoints(i,1)),round(Zpoints(i,3)))=100;
% end
% Zstack_filt2=imgaussfilt3(Zstack_filt2,5);
%
%
% confocal_filt2=confocal_filt2.*confocal_filt;
% Zstack_filt2=Zstack_filt2.*Zstack_filt;

% pixel information extraction from each cell
figure
patchradius=[150 150 15];
for i=1:size(Spoints,1)
    newpatchradius=patch_bound(Spoints(i,[2 1 3]),patchradius(1:3),size(confocal_filt));
    S_patch{i}=max(confocal_filt(Spoints(i,2)-newpatchradius(1):Spoints(i,2)+newpatchradius(1),Spoints(i,1)-newpatchradius(2):Spoints(i,1)+newpatchradius(2),Spoints(i,3)-newpatchradius(3):Spoints(i,3)+newpatchradius(3)),[],3);
    S_patchsize{i}=size(S_patch{i});
    S_patch{i}=padarray(S_patch{i},floor([patchradius(1)+0.5 patchradius(2)+0.5]-floor(size(S_patch{i})/2)),'both');
    if size(S_patch{i},1)==2*patchradius(1)
        S_patch{i}(end+1,:)=0;
    end
        if size(S_patch{i},2)==2*patchradius(2)
        S_patch{i}(:,end+1)=0;
    end
    imagesc(S_patch{i});axis equal;axis off
    drawnow
end
sIM=[];
t=0;
for i=1:round(sqrt(length(S_patch)))
    tmp=[];
    for j=1:ceil(sqrt(length(S_patch)))
        t=t+1;
        try;tmp=[tmp linhistmatch(S_patch{t}(140:160,140:160),S_patch{1}(140:160,140:160),50,'regular')];
        catch
            tmp=[tmp zeros(21,21)];
        end
    end
    sIM=[sIM;tmp];
end

for i=1:length(S_patch)
    for j=1:length(S_patch)
        Spatchsim(i,j)=corr(vec(S_patch{i}(140:160,140:160)),vec(S_patch{j}(140:160,140:160)));
    end
end

for i=1:size(Zpoints,1)
    newpatchradius=patch_bound(Zpoints(i,[2 1 3]),patchradius(1:3),size(Zstack_filt));
    Z_patch{i}=max(Zstack_filt(Zpoints(i,2)-newpatchradius(1):Zpoints(i,2)+newpatchradius(1),Zpoints(i,1)-newpatchradius(2):Zpoints(i,1)+newpatchradius(2),Zpoints(i,3)-newpatchradius(3):Zpoints(i,3)+newpatchradius(3)),[],3);
    Z_patchsize{i}=size(Z_patch{i});
    Z_patch{i}=padarray(Z_patch{i},floor([patchradius(1)+0.5 patchradius(2)+0.5]-floor(size(Z_patch{i})/2)),'both');
    if size(Z_patch{i},1)==2*patchradius(1)
        Z_patch{i}(end+1,:)=0;
    end
        if size(Z_patch{i},2)==2*patchradius(2)
        Z_patch{i}(:,end+1)=0;
    end
    imagesc(Z_patch{i});axis equal;axis off
    drawnow
end

zIM=[];
t=0;
for i=1:round(sqrt(length(Z_patch)))
    tmp=[];
    for j=1:ceil(sqrt(length(Z_patch)))
        t=t+1;
        try;tmp=[tmp linhistmatch(Z_patch{t}(140:160,140:160),S_patch{1}(140:160,140:160),50,'regular')];
        catch
            tmp=[tmp zeros(21,21)];
        end
    end
    zIM=[zIM;tmp];
end


for i=1:length(Z_patch)
    for j=1:length(Z_patch)
        Zpatchsim(i,j)=corr(vec(Z_patch{i}(140:160,140:160)),vec(Z_patch{j}(140:160,140:160)));
    end
end
%
% angles=linspace(0,360,20);
% bestsofar=0;
% for i=1:length(S_patch)
%
%     for j=1:length(Z_patch)
%         for a=1:length(angles)
%             [S_patch_trans,D_im(i,j),~]=linhistmatch(S_patch{i},Z_patch{j},200,'regular');
%             Z_rot=log1p(imrotate(Z_patch{j},angles(a),'bilinear','crop'));
%             XC(i,j,a)=vec(log1p(S_patch_trans))'*vec(Z_rot);
%             if XC(i,j,a)>bestsofar;
%                 bestsofar=XC(i,j,a);
%                 subplot(1,3,1)
%                 imagesc(S_patch{i});
%                 subplot(1,3,2)
%                 imagesc(Z_rot);
%                 subplot(1,3,3)
%                 imagesc(S_patch{i}.*Z_rot);
%                 drawnow
%             end
%             [i j a]
%         end
%     end
% end
% [vol_dist]=knn_distance(Spoints,Zpoints,10,30);

%% Point set registration
figure(3)
subplot(1,2,1)
hold on
plot3(Spoints(:,1),Spoints(:,2),Spoints(:,3),'b.','MarkerSize',20);axis equal;axis tight
plot3(Zpoints(:,1),Zpoints(:,2),Zpoints(:,3),'r.','MarkerSize',20);axis equal;axis tight
legend('Confocal','Zstack');
set(gca,'FontWeight','bold','FontSize',20,'TickLength',[0 0]);set(gcf,'Color','w');


%% homage to all the previous methods that were garbage or slow (lol)
% [bhat,P,inlier_set,q]=rrwoc_wahba(Spoints,Zpoints,0.9,30,1,2000);

% [vol_dist,P]=volumetric_distance_2(Spoints,Zpoints,10);

% [vol_dist]=knn_distance(Spoints,Zpoints,50,30);
% M=munkres(-vol_dist);
%
% [R,T,beta]=weighted_wahba(Spoints,M*Zpoints,max(M.*exp(vol_dist./2*0.1^2),[],2));
%

% [pl_dist]=procrustes_levenstein_distance(Spoints,Zpoints,radius_parameter,tightness_paramter,max_depth_parameter);
% [s_idx,z_idx]=find(pl_dist==min(pl_dist(:)));
% [R,T]=procrustes_levenstein_distance_decoder(Spoints,Zpoints,s_idx,z_idx,radius_parameter,tightness_paramter,max_depth_parameter);

[pl_dist,maximal_rotation_set,satisfies_anglelimit,maximal_rotation,maximal_rotation_founders]=maximal_rotation_group(Spoints,Zpoints,radius_parameter,tightness_parameter,max_depth_parameter,anglelimit);
[s_idx,z_idx]=find((pl_dist.*satisfies_anglelimit)==max(vec(pl_dist.*satisfies_anglelimit)));
clear R;for t=1:size(s_idx,1);R{t}=maximal_rotation{s_idx(t),z_idx(t)};end

%% Visualizing registration matches to certify by eye
for t=1:length(R)
    disp(['Showing solution ' num2str(t) ' out of ' num2str(length(R)) ' - number of matches: ' num2str(pl_dist(s_idx(t),z_idx(t)))]);
    bhat=R{t};
    
    
    
    Zpointshat=[Spoints ones(size(Spoints,1),1)]*bhat;
    
%     [P.matched_cells_y,P.matched_cells_x]=find(pdist2(Zpoints,Zpointshat)<= 2*tightness_parameter);
    P.matched_cells_y=maximal_rotation_set{s_idx(t),z_idx(t)}(:,2);
    P.matched_cells_x=maximal_rotation_set{s_idx(t),z_idx(t)}(:,1);
    figure(3)
    subplot(1,2,2)
    hold on
    plot3(Zpointshat(:,1),Zpointshat(:,2),Zpointshat(:,3),'b.','MarkerSize',20);axis equal;axis tight
    plot3(Zpoints(:,1),Zpoints(:,2),Zpoints(:,3),'r.','MarkerSize',20);axis equal;axis tight
    legend('Registered Confocal','Zstack');
    set(gca,'FontWeight','bold','FontSize',20,'TickLength',[0 0]);set(gcf,'Color','w');
    
    
    confocal_histmatched=linhistmatch(log1p(confocal_filt),minmax(max(log1p(Zstack_filt),[],3)),200,'regular');
    tic
    tform = affine3d([[bhat(1:3,1:3);bhat(4,:)] [0;0;0;1]]);
    confocal_moved=imwarp(confocal_histmatched,tform,'Outputview',imref3d(size(Zstack_filt)));
    toc
    
    
    figure(4)
    
    A=max(confocal_moved(:,:,max(Zpoints(z_idx(t),3)-5,1):min(Zpoints(z_idx(t),3)+5,size(confocal_moved,3))),[],3);
    B=max(minmax(log1p(Zstack_filt(:,:,max(Zpoints(z_idx(t),3)-5,1):min(Zpoints(z_idx(t),3)+5,size(confocal_moved,3))))),[],3);
    
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
    title('Confocal/Multiplex/Z-stack (left/middle/right)');
    
    figure(5)
    
    newpatchradius=patch_bound(Zpoints(z_idx(1),:),[200 200 20],size(Zstack_filt));
    Ap=max(confocal_moved(Zpoints(z_idx(1),2)-newpatchradius(2):Zpoints(z_idx(1),2)+newpatchradius(2),Zpoints(z_idx(1),1)-newpatchradius(1):Zpoints(z_idx(1),1)+newpatchradius(1),max(Zpoints(z_idx(t),3)-5,1):min(Zpoints(z_idx(t),3)+5,size(confocal_moved,3))),[],3);
    Bp=max(minmax(log1p(Zstack_filt(Zpoints(z_idx(1),2)-newpatchradius(2):Zpoints(z_idx(1),2)+newpatchradius(2),Zpoints(z_idx(1),1)-newpatchradius(1):Zpoints(z_idx(1),1)+newpatchradius(1),max(Zpoints(z_idx(t),3)-5,1):min(Zpoints(z_idx(t),3)+5,size(confocal_moved,3))))),[],3);
    
    Ip(:,:,1)=[Ap Ap.*Bp zeros(size(Bp))];
    Ip(:,:,2)=[zeros(size(Ap)) zeros(size(Ap)) zeros(size(Bp))];
    Ip(:,:,3)=[zeros(size(Ap)) Ap.*Bp Bp];
    imagesc(Ip)
    title('Zoom of Confocal/Multiplex/Z-stack (left/middle/right)');
    pause
end
