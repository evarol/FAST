function [vol_dist,dx,dy,tx,ty,vx,vy]=volumetric_distance(X,Y);

[DX,xnn]=sort(squareform(pdist(X)),2,'ascend');
[DY,ynn]=sort(squareform(pdist(Y)),2,'ascend');
for i=1:size(X,1)
    dx(i,:)=DX(i,2:4);
end
for j=1:size(Y,1)
    dy(j,:)=DY(j,2:4);
end
for i=1:size(X,1)
    tx(i,1)=heron(X(xnn(i,[1 2 3]),:));
    tx(i,2)=heron(X(xnn(i,[1 2 4]),:));
    tx(i,3)=heron(X(xnn(i,[1 3 4]),:));
end
tx=sort(tx,2,'ascend');
for j=1:size(Y,1)
    ty(j,1)=heron(Y(ynn(j,[1 2 3]),:));
    ty(j,2)=heron(Y(ynn(j,[1 2 4]),:));
    ty(j,3)=heron(Y(ynn(j,[1 3 4]),:));
end
ty=sort(ty,2,'ascend');

for i=1:size(X,1)
    [~,vx(i,1)]=convhulln(X(xnn(i,1:4),:));
end
for j=1:size(Y,1)
    [~,vy(j,1)]=convhulln(Y(ynn(j,1:4),:));
end


vol_dist(:,:,1)=pdist2(dx,dy);
vol_dist(:,:,2)=pdist2(tx,ty);
vol_dist(:,:,3)=pdist2(vx,vy);
end


function area=heron(X)
a=norm(X(1,:)-X(2,:));
b=norm(X(1,:)-X(3,:));
c=norm(X(2,:)-X(3,:));
s = (a+ b + c)/2;
area = sqrt(s * (s-a) * (s-b) * (s-c));
end