clear all
clc
close all
vec=@(x)(x(:));
disp('Select confocal point clouds(s)...');
[confocalfile,confocalpath] = uigetfile('*','Select files to load','MultiSelect','on');
disp('Select zstack point cloud...');
[zstackfile,zstackpath] = uigetfile('*','Select files to load','MultiSelect','on');


if ~isequal(confocalfile,0)
    if ~iscell(confocalfile)
        tmp{1}=confocalfile;
        confocalfile=tmp;
    end
    for i=1:length(confocalfile)
        tmp2=load([confocalpath confocalfile{i}]);
        Spoints{i}=tmp2.Spoints;
    end
end


if ~isequal(zstackfile,0)
    if ~iscell(zstackfile)
        tmp{1}=zstackfile;
        zstackfile=tmp;
    end
    
    for i=1:length(zstackfile)
        tmp2=load([zstackpath zstackfile{i}]);
        Zpoints{i}=tmp2.Zpoints;
    end
end

maxdist=200;
indelcost=10;
zmax=30;
anglelimit=20;
for i=1:length(Spoints)
    for j=1:length(Zpoints)
        X=Spoints{i};
        Y=Zpoints{j};
        [max_rotation_dist{i,j},maximal_rotation_set{i,j},satisfies_anglelimit{i,j},maximal_rotation{i,j},maximal_rotation_founders{i,j}]=maximal_rotation_group(X,Y,maxdist,indelcost,zmax,anglelimit);
        bestreg(i,j)=max(vec(max_rotation_dist{i,j}.*satisfies_anglelimit{i,j}));
    end
end
