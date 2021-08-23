function view_annotations

minmax = @(x)((x-min(x(:)))./max(x(:)-min(x(:))));
addpath ./tiff_loading/jsonlab/
addpath ./tiff_loading/utilities/
addpath(genpath('./tiff_loading/Fiji.app'));
javaaddpath('./tiff_loading/Fiji.app/mij.jar')
disp('Select image to annotate...');
[file,path] = uigetfile('*','Select image file to load','MultiSelect','on');
imfull=double(load_tiff([path file]));
if size(imfull,4)>1
    imfull=imfull(:,:,:,2);
end
[file,path] = uigetfile('*','Select annotation file to load','MultiSelect','on');
load([path file]);


close all

imagesc(max(log1p(imfull),[],3));
hold on
plot(points(:,1),points(:,2),'r.','MarkerSize',20);
title('Final top down view of annotations');

composite_image=repmat(minmax(log1p(imfull)),[1 1 1 3]);
for i=1:size(points,1)
    composite_image(points(i,2)-5:points(i,2)+5,points(i,1)-5:points(i,1)+5,points(i,3),1)=1;
end
figure
imshow3d(composite_image);
title('Final annotations in 3D-scroller');

end