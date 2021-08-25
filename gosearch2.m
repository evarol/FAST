function [max_pl_dist]=gosearch2(X,Y,maxdist,indelcost,zmax,anglelimit,gridsize);
anglelimit=deg2rad([0 360;-anglelimit anglelimit;-anglelimit anglelimit]);
vec=@(x)(x(:));
DX=squareform(pdist(X));
DY=squareform(pdist(Y));
for i=1:size(X,1)
    for j=1:size(X,1)
        if abs(X(i,3)-X(j,3))>zmax
            DX(i,j)=Inf;
        end
    end
end
for i=1:size(Y,1)
    for j=1:size(Y,1)
        if abs(Y(i,3)-Y(j,3))>zmax
            DY(i,j)=Inf;
        end
    end
end

[SX,nnx]=sort(DX,2,'ascend');
[SY,nny]=sort(DY,2,'ascend');
kx=sum(DX<=maxdist,2);
ky=sum(DY<=maxdist,2);
pl_dist=zeros(size(X,1),size(Y,1));
bestsofar=0;
for i=1:size(X,1)
    for j=1:size(Y,1)
        [i j]
        [~,~,~,xhat,~]=nw(SX(i,1:kx(i))',SY(j,1:ky(j))',indelcost);
        max_pl_dist(i,j)=sum(~isnan(xhat));
    end
end




