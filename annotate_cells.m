function annotate_cells(downsampling_factor)
if nargin<1
    downsampling_factor=1;
end
minmax = @(x)((x-min(x(:)))./max(x(:)-min(x(:))));
addpath ./tiff_loading/jsonlab/
addpath ./tiff_loading/utilities/
addpath(genpath('./tiff_loading/Fiji.app'));
javaaddpath('./tiff_loading/Fiji.app/mij.jar')
disp('Select image to annotate...');
[file,path] = uigetfile('*','Select file to load','MultiSelect','on');
im=double(load_tiff([path file]));

if size(im,4)>1
    im=im(:,:,:,2);
end


current_slice = 1;
imfull=im;
im=imresize3(log1p(im),[size(im,1) size(im,2) size(im,3)/downsampling_factor]);
keep_going=1;
points=[];
figure('units','normalized','outerposition',[0 0 0.5 1])
while keep_going==1
    cla;imagesc(im(:,:,current_slice));
    if ~isempty(points);hold on;plot(points(:,1),points(:,2),'m.','MarkerSize',20);end
    text(40,40,[num2str(current_slice) '/' num2str(size(im,3))],'FontSize',30,'FontWeight','bold','Color','w');
    title('Left click to select cells, right click to end annotation, delete to erase previous annotation');
    [x,y]=getpts;
    points=[points;[x(1:end-1) y(1:end-1) repmat(current_slice*downsampling_factor,[size(x,1)-1 1])]];
    title('Next slide? (Enter), Re-annotate? (R+Enter), Go to previous slice? (P+Enter), Finished? (F+Enter)');
    prompt=input('Next slide? (Enter), Re-annotate? (R+Enter), Go to previous slice? (P+Enter), Finished? (F+Enter)','s');
    if isempty(prompt)
        current_slice=current_slice+1;
        if current_slice>size(im,3)
            keep_going=0;
        end
    end
    if strcmpi(prompt,'P')
        points(points(:,3)==current_slice,:)=[];
        current_slice=max(current_slice-1,1);
        points(points(:,3)==current_slice,:)=[];
    end
    
    if strcmpi(prompt,'R')
        points(points(:,3)==current_slice,:)=[];
    end
    if strcmpi(prompt,'F')
        keep_going=0;
    end
end

if isempty(points)
    return
end
close all

imagesc(max(log1p(imfull),[],3));
hold on
plot(points(:,1),points(:,2),'r.','MarkerSize',20);
title('Final top down view of annotations');

composite_image=repmat(minmax(log1p(imfull)),[1 1 1 3]);
for i=1:size(points,1)
    composite_image(points(i,2)-5:points(i,2)+5,points(i,1)-5:points(i,1)+5,max(points(i,3)-downsampling_factor/2,1):min(points(i,3)+downsampling_factor/2,size(imfull,3)),1)=1;
end
figure
imshow3d(composite_image);
title('Final annotations in 3D-scroller');
names=strsplit(file,'.');
save([path names{1} '_annotations.mat'],'points');

end