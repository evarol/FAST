function [pl_dist,ransac]=procrustes_levenstein_ransac(X,Y,maxdist,indelcost,zmax)
%% finds the levenstein distance between the k-nearest neighbor distances of two sets of points
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
pl_dist=inf(size(X,1),size(Y,1));
ransac=zeros(size(X,1),size(Y,1),2);
for i=1:size(X,1)
    for j=1:size(Y,1)
        [i j]
        if or(kx(i)<4,ky(j)<4)
            pl_dist(i,j)=Inf;
        else
            x_candidates=combnk(2:kx(i),3);
            for ii=1:size(x_candidates,1)
                [T,~,~,xhat,match]=nw(SX(i,[1 x_candidates(ii,:)])',SY(j,1:ky(j))',indelcost);
                if any(isnan(xhat))
                    pl_dist(i,j)=Inf;
                else
                    Xmatch=[X(i,:);X(nnx(i,x_candidates(ii,:)),:)];
                    Ymatch=Y(nny(j,flipud(match(:,2))),:);
                    [R,T,~]=wahba(Xmatch,Ymatch);
                    angles=rad2deg(rotm2eul(R));
                    if and(abs(angles(2))<30,abs(angles(3))<30)
                       tmp=(pdist2(X*R+T,Y)<=indelcost);
                       ransac(:,:,1)=ransac(:,:,1)+tmp;
                       ransac(:,:,2)=max(ransac(:,:,2),tmp*sum(vec(tmp)));
                    pl_dist(i,j)=min(pl_dist(i,j),-sum(vec(tmp)));
                    end
                end
            end
        end
    end   
end