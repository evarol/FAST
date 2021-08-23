function [Ysubset,hardGW]=volumetric_sampling(Xsubset,Yleft,Yfirst)



%% Sampling 1st point -- satisfy null volume -- anything works here
if nargin<3
    Y=Yleft;
r=randperm(size(Y,1));
Ysubset(1,:)=Y(r(1),:);
Yleft=Y(r(2:end),:);
else
    if length(Yfirst)==1
    Ysubset(1,:)=Yleft(Yfirst,:);
    Yleft=Yleft((1:size(Yleft,1))~=Yfirst,:);
    else
        Ysubset(1,:)=Yfirst;
    end
end

%% Sampling 2nd point --- 2nd point must yield a distance to the first point that is included in the pairwise distances in Xsubset

Dcandidate=pdist2(Yleft,Ysubset);

Dactual = pdist(Xsubset)';

Dcomp = min(pdist2(Dcandidate,Dactual),[],2);

[~,idx]=min(Dcomp);

Ysubset(2,:)=Yleft(idx,:);
Yleft(idx,:)=[];

%% Sampling 3rd point --- 3rd point must yield an area to the first 2 points that is included in the areas in Xsubset

Aactual(1,1)=heron(Xsubset([1 2 3],:));
Aactual(2,1)=heron(Xsubset([1 2 4],:));
Aactual(3,1)=heron(Xsubset([4 2 3],:));


for i=1:size(Yleft,1)
    Acandidate(i,1)=heron([Ysubset;Yleft(i,:)]);
end

Acomp = min(pdist2(Acandidate,Aactual),[],2);

[~,idx]=min(Acomp);

Ysubset(3,:)=Yleft(idx,:);
Yleft(idx,:)=[];


%% Sampling 4th point --- 4th point must yield a volume to the first 3 points that is similar to the volume of Xsubset

try;[~,Vactual]=convhull(Xsubset);
catch
    Vactual=0;
end

for i=1:size(Yleft,1)
try;    [~,Vcandidate(i,1)]=convhull([Ysubset;Yleft(i,:)]);
catch
    Vcandidate(i,1)=0;
end
end

Vcomp = min(pdist2(Vcandidate,Vactual),[],2);
[~,idx]=min(Vcomp);
Ysubset(4,:)=Yleft(idx,:);
Yleft(idx,:)=[];


%% Do combinatorial graph matching - 24 combos - just do it

P=perms([1 2 3 4]);

for i=1:size(P,1)
    Pcandidate(i,1)=norm(pdist(Xsubset)-pdist(Ysubset(P(i,:),:)));
end
[hardGW,idx]=min(Pcandidate);

phat = P(idx,:);
Ysubset=Ysubset(phat,:);



end

function area=heron(X)
a=norm(X(1,:)-X(2,:));
b=norm(X(1,:)-X(3,:));
c=norm(X(2,:)-X(3,:));
s = (a+ b + c)/2;
area = sqrt(s * (s-a) * (s-b) * (s-c));
end