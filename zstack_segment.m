function [Zpoints,Zstats,Zstack_filt]=zstack_segment(zstack,z,sigma)
%% auxillary functions
minmax = @(x)((x-min(x(:)))./max(x(:)-min(x(:))));
minmaxnnz = @(x)(max((x-min(x(x>0))+1)./max(x(:)-min(x(x>0))+1),0));
vec=@(x)(x(:));


Izstack=log1p(zstack);
Bzstack=zeros(size(Izstack));
Zstack_filt=zeros(size(zstack));
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