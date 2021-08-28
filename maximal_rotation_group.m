function [max_rotation_dist,maximal_rotation_set,satisfies_anglelimit,maximal_rotation,maximal_rotation_founders]=maximal_rotation_group(X,Y,maxdist,indelcost,zmax,anglelimit)
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
max_rotation_dist=zeros(size(X,1),size(Y,1));
maximal_rotation_set=cell(size(X,1),size(Y,1));
maximal_rotation_founders=cell(size(X,1),size(Y,1));
satisfies_anglelimit=zeros(size(X,1),size(Y,1));
maximal_rotation=cell(size(X,1),size(Y,1));
numsearch=size(X,1)*size(Y,1);t=0;tic;
for i=1:size(X,1)
    for j=1:size(Y,1)
        t=t+1;
        if and(kx(i)>2,ky(i)>2)
            D=(pdist2(SX(i,1:kx(i))',SY(j,1:ky(j))')<=indelcost);
            xmatch=find(any(D==1,2));
            if length(xmatch)>2
                x_candidates=nchoosek(xmatch(2:end),2);
                for ii=1:size(x_candidates,1)
                    x_set=[i nnx(i,x_candidates(ii,1)) nnx(i,x_candidates(ii,2))]';
                    y_candidates1=find(D(x_candidates(ii,1),:)==1);
                    y_candidates2=find(D(x_candidates(ii,2),:)==1);
                    for j1=1:length(y_candidates1)
                        for j2=1:length(y_candidates2)
                            y_set=[j nny(j,y_candidates1(j1)) nny(j,y_candidates2(j2))]';
                            if isrotation(DX(x_set,x_set),DY(y_set,y_set),indelcost)
                                [x_new,y_new]=growrotation(DX,DY,x_set,y_set,indelcost,X,Y,anglelimit);
                                if length(x_new)>max_rotation_dist(i,j)
                                    max_rotation_dist(i,j)=length(x_new);
                                    maximal_rotation_set{i,j}=[x_new,y_new];
                                    maximal_rotation_founders{i,j}=[x_set,y_set];
                                    [angles,R,T]=check_angle(X(x_new,:),Y(y_new,:));
                                    maximal_rotation{i,j}=[R;T];
                                    satisfies_anglelimit(i,j)=and(abs(angles(2))<anglelimit,abs(angles(3))<anglelimit);
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    clc
    fprintf(['Searching for transformations (' num2str(t) '/' num2str(numsearch) ')...\n']);
    fprintf(['\n' repmat('.',1,50) '\n\n'])
    for tt=1:round(t*50/(numsearch))
        fprintf('\b|\n');
    end
    TT=toc;
    disp(['Time elapsed (minutes): ' num2str(TT/60) ' Time remaining (minutes): ' num2str((numsearch-t)*(TT/t)*(1/60))]);
end
end

function out=isrotation(DX,DY,indelcost)

if max(abs(DX(:)-DY(:)))<indelcost
    out=1;
else
    out=0;
end
end

function [x_set,y_set]=growrotation(DX,DY,x_founders,y_founders,indelcost,X,Y,anglelimit);
vec=@(x)(x(:));
x_candidates=setdiff((1:size(DX,1))',x_founders);
y_candidates=setdiff((1:size(DY,1))',y_founders);
DXDY=pdist2(DX(x_candidates,x_founders),DY(y_candidates,y_founders));

[x_sel,y_sel]=find(DXDY<=indelcost);
max_set=0;
max_closeness=Inf;
x_fin_sel=[];
y_fin_sel=[];

%% its choosing more stuff - going away from founders
for i=1:size(x_sel,1)
    R=wahba(X([x_founders;x_candidates(x_sel(i))],:),Y([y_founders;y_candidates(y_sel(i))],:));
    angles=rad2deg(rotm2eul(R));
    if and(abs(angles(2))<anglelimit,abs(angles(3))<anglelimit)
        subDXDY=pdist2(DX(x_candidates,[x_founders;x_candidates(x_sel(i))]),DY(y_candidates,[y_founders;y_candidates(y_sel(i))]));
        [x_subsel,y_subsel]=find(subDXDY<=indelcost);
        set_size= sum(and(isunique(x_subsel),isunique(y_subsel)));
        set_closeness= max(vec(subDXDY(x_subsel,y_subsel)));
        if set_size>=max_set
            max_set=set_size;
            if set_closeness<max_closeness
                max_closeness=set_closeness;
                x_fin_sel=[x_candidates(x_subsel)];
                y_fin_sel=[y_candidates(y_subsel)];
            end
            
        end
    end
end
x_set=[x_founders;x_fin_sel];
y_set=[y_founders;y_fin_sel];
end

function u=isunique(x)
[x,idx]=sort(x);
tmp=ones(size(x,1),1);
u=ones(size(x,1),1);
for i=1:size(x,1)-1
    if x(i,1)==x(i+1,1)
        tmp(i)=0;
        tmp(i+1)=0;
    end
end
u(idx)=tmp;
end

function [angles,R,T]=check_angle(X,Y)
[R,T]=wahba(X,Y);
angles=rad2deg(rotm2eul(R));
end