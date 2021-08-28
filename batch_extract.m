clear all
clc
close all

%% dependencies to load tiffs/nd2s
addpath ./tiff_loading/jsonlab/
addpath ./tiff_loading/utilities/
addpath(genpath('./tiff_loading/Fiji.app'));
javaaddpath('./tiff_loading/Fiji.app/mij.jar')


%% auxillary functions
minmax = @(x)((x-min(x(:)))./max(x(:)-min(x(:))));
minmaxnnz = @(x)(max((x-min(x(x>0))+1)./max(x(:)-min(x(x>0))+1),0));
vec=@(x)(x(:));

%% image parameters
confocal_scale = [1272.7922/2048 1272.7922/2048 44/22]; %x,y,z pixel dimension - user specified. VERY IMPORTANT to have real spatial dimensions
zstack_scale = [1000/700 1000/700 2]; %x,y,z pixel dimension - user specified. VERY IMPORTANT to have real spatial dimensions

%% Detection parameters
z=6; %z-score to distinguish bright foreground from background
sigma=2; %smoothing kernel to fill in gaps (to avoid spotty cells)

%% load data
disp('Select confocal image(s)...');
[confocalfile,confocalpath] = uigetfile('*','Select files to load','MultiSelect','on');
disp('Select zstack image(s)...');
[zstackfile,zstackpath] = uigetfile('*','Select files to load','MultiSelect','on');


if ~isequal(confocalfile,0)
    if ~iscell(confocalfile)
        tmp{1}=confocalfile;
        confocalfile=tmp;
    end
    
    t=0;numsearch=length(confocalfile);tic;
    for i=1:length(confocalfile)
        close all
        t=t+1;
        confocal=double(load_tiff([confocalpath confocalfile{i}]));
        confocal=double(confocal(:,:,:,2));
        confocal_size = size(confocal).*confocal_scale;
        confocal=max(imresize3(confocal,confocal_size),0);
        [Spoints,Sstats,confocal_filt]=confocal_segment(confocal,z,sigma);
        names=strsplit(confocalfile{i},'.');
        save([confocalpath names{1} '_point_set.mat'],'Spoints');
        clc
        fprintf(['Extracting confocal point sets (' num2str(t) '/' num2str(numsearch) ')...\n']);
        fprintf(['\n' repmat('.',1,50) '\n\n'])
        for tt=1:round(t*50/(numsearch))
            fprintf('\b|\n');
        end
        TT=toc;
        disp(['Time elapsed (minutes): ' num2str(TT/60) ' Time remaining (minutes): ' num2str((numsearch-t)*(TT/t)*(1/60))]);
        MIJ.run("Clear Results");
        MIJ.closeAllWindows
    end
end


if ~isequal(zstackfile,0)
    if ~iscell(zstackfile)
        tmp{1}=zstackfile;
        zstackfile=tmp;
    end
    
    t=0;numsearch=length(zstackfile);tic;
    for i=1:length(zstackfile)
        close all
        t=t+1;
        zstack=double(load_tiff([zstackpath zstackfile{i}]));
        zstack_size = size(zstack).*zstack_scale;
        zstack=max(imresize3(zstack,zstack_size),0);
        [Zpoints,Zstats,Zstack_filt]=zstack_segment(zstack,z,sigma);
        names=strsplit(zstackfile{i},'.');
        save([zstackpath names{1} '_point_set.mat'],'Zpoints');
        clc
        fprintf(['Extracting confocal point sets (' num2str(t) '/' num2str(numsearch) ')...\n']);
        fprintf(['\n' repmat('.',1,50) '\n\n'])
        for tt=1:round(t*50/(numsearch))
            fprintf('\b|\n');
        end
        TT=toc;
        disp(['Time elapsed (minutes): ' num2str(TT/60) ' Time remaining (minutes): ' num2str((numsearch-t)*(TT/t)*(1/60))]);
        MIJ.run("Clear Results");
        MIJ.closeAllWindows
    end
end
