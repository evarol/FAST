function [pl_dist,R,T]=gosearch(X,Y,sigma,anglelimit,gridsize);
anglelimit=[0 360;-anglelimit anglelimit;-anglelimit anglelimit];
pl_dist=zeros(size(X,1),size(Y,1));
bestsofar=0;
tic;t=0;
numsearch=size(X,1)*size(Y,1);
for i=1:size(X,1)
    for j=1:size(Y,1)
        t=t+1;
        [bestmatch,Rbest]=brute_rot_3dxcorr(X+Y(j,:)-X(i,:),Y,sigma,anglelimit,gridsize);
        pl_dist(i,j)=bestmatch;
        if bestmatch>bestsofar
            bestsofar=bestmatch;
            R=Rbest;
            T=Y(j,:)-X(i,:);
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
end

