function [Rout,Tout]=procrustes_levenstein_distance_decoder(X,Y,x_idx,y_idx,maxdist,indelcost,zmax)
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

t=0;
for iii=1:length(x_idx)
    i=x_idx(iii);
    j=y_idx(iii);
    [i j]
    bestsofar=Inf;
    t=t+1;
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
                [R,T,beta]=wahba(Xmatch,Ymatch);
                pl_dist(i,j)=min(pl_dist(i,j),-sum(vec(pdist2(X*R+T,Y)<=indelcost)));
                if pl_dist(i,j)<bestsofar
                    bestsofar=pl_dist(i,j);
                    Rout{t}=R;Tout{t}=T;
                end
            end
        end
    end
end