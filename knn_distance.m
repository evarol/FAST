function [vol_dist]=knn_distance(X,Y,k,indelcost);
%% finds the levenstein distance between the k-nearest neighbor distances of two sets of points
DX=sort(squareform(pdist(X)),2,'ascend');
DY=sort(squareform(pdist(Y)),2,'ascend');
vol_dist=zeros(size(X,1),size(Y,1));
for i=1:size(X,1)
    for j=1:size(Y,1)
       [T,P,skip,xhat]=nw(DX(i,1:k)',DY(j,1:k)',indelcost);
       vol_dist(i,j)=T(end-1,end-1);
       [i j]
    end
end


end