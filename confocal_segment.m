function [Spoints,Sstats,confocal_filt]=confocal_segment(confocal,z,sigma)
%% auxillary functions
minmax = @(x)((x-min(x(:)))./max(x(:)-min(x(:))));
minmaxnnz = @(x)(max((x-min(x(x>0))+1)./max(x(:)-min(x(x>0))+1),0));
vec=@(x)(x(:));

Iconfocal=log1p(confocal);
Bconfocal=zeros(size(Iconfocal));
confocal_filt=zeros(size(confocal));
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