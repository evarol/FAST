function [lev_dist,num_match,proc_distance,hull_distance,kx,ky]=max_dist_nn_distance(X,Y,maxdist,indelcost,visualization)
%% finds the levenstein distance between the k-nearest neighbor distances of two sets of points
DX=squareform(pdist(X));
DY=squareform(pdist(Y));
for i=1:size(X,1)
    for j=1:size(X,1)
        if abs(X(i,3)-X(j,3))>15
            DX(i,j)=Inf;
        end
    end
end
for i=1:size(Y,1)
    for j=1:size(Y,1)
        if abs(Y(i,3)-Y(j,3))>15
            DY(i,j)=Inf;
        end
    end
end

[SX,nnx]=sort(DX,2,'ascend');
[SY,nny]=sort(DY,2,'ascend');
kx=sum(DX<=maxdist,2);
ky=sum(DY<=maxdist,2);
lev_dist=zeros(size(X,1),size(Y,1));
num_match=zeros(size(X,1),size(Y,1));
proc_distance=zeros(size(X,1),size(Y,1));
hull_distance=zeros(size(X,1),size(Y,1));
for i=1:size(X,1)
    for j=1:size(Y,1)
        [T,~,~,xhat,match]=nw(SX(i,1:kx(i))',SY(j,1:ky(j))',indelcost);
        [i j]
        lev_dist(i,j)=T(end,end);
        num_match(i,j)=(sum(~isnan(xhat)));
        tmp=X(nnx(i,1:kx(i)),:);
        Xmatch=tmp(match(:,1),:);
        tmp=Y(nny(j,1:ky(j)),:);
        Ymatch=tmp(match(:,2),:);
        [R,T,~]=wahba(Xmatch,Ymatch);
        proc_distance(i,j)=mean(norms(Xmatch*R + T - Ymatch,2));
        try;[~,vx]=convhulln(Xmatch);[~,vy]=convhulln(Ymatch);
        hull_distance(i,j)=abs(vx-vy);catch; hull_distance(i,j)=Inf;end
       if visualization==1
            [SX(i,match(:,1))' SY(j,match(:,2))']
            if any(abs(SX(i,match(:,1))'-SY(j,match(:,2))')>indelcost)
                warning('Violation');
            end
            subplot(1,3,1)
            cla
            hold on
            plot3(X(nnx(i,1:kx(i)),1),X(nnx(i,1:kx(i)),2),X(nnx(i,1:kx(i)),3),'b.','MarkerSize',20);
            for k=1:kx(i)
                plot3([X(nnx(i,1),1) X(nnx(i,k),1)],[X(nnx(i,1),2) X(nnx(i,k),2)],[X(nnx(i,1),3) X(nnx(i,k),3)],'k-','LineWidth',2);
                text(X(nnx(i,k),1),X(nnx(i,k),2),X(nnx(i,k),3),num2str(SX(i,k)),'Color','k','FontWeight','bold','FontSize',20);
            end
            axis equal;axis tight;
            title('X-set');
            
            subplot(1,3,2)
            cla
            hold on
            plot3(Y(nny(j,1:ky(j)),1),Y(nny(j,1:ky(j)),2),Y(nny(j,1:ky(j)),3),'r.','MarkerSize',20);
            for k=1:ky(j)
                plot3([Y(nny(j,1),1) Y(nny(j,k),1)],[Y(nny(j,1),2) Y(nny(j,k),2)],[Y(nny(j,1),3) Y(nny(j,k),3)],'k-','LineWidth',2);
                text(Y(nny(j,k),1),Y(nny(j,k),2),Y(nny(j,k),3),num2str(SY(j,k)),'Color','k','FontWeight','bold','FontSize',20);
            end
            axis equal;axis tight;title('Y-set');
            
            subplot(1,3,3)
            cla
            hold on
            plot3(X(nnx(i,1:kx(i)),1),X(nnx(i,1:kx(i)),2),X(nnx(i,1:kx(i)),3),'b.','MarkerSize',20);
            for k=1:kx(i)
                if ismember(k,match(:,1))
                    plot3([X(nnx(i,1),1) X(nnx(i,k),1)],[X(nnx(i,1),2) X(nnx(i,k),2)],[X(nnx(i,1),3) X(nnx(i,k),3)],'k-','LineWidth',2);
                    text(X(nnx(i,k),1),X(nnx(i,k),2),X(nnx(i,k),3),num2str(SX(i,k)),'Color','k','FontWeight','bold','FontSize',20);
                end
            end
            
            
            plot3(Y(nny(j,1:ky(j)),1),Y(nny(j,1:ky(j)),2),Y(nny(j,1:ky(j)),3),'r.','MarkerSize',20);
            for k=1:ky(j)
                if ismember(k,match(:,2))
                    plot3([Y(nny(j,1),1) Y(nny(j,k),1)],[Y(nny(j,1),2) Y(nny(j,k),2)],[Y(nny(j,1),3) Y(nny(j,k),3)],'k-','LineWidth',2);
                    text(Y(nny(j,k),1),Y(nny(j,k),2),Y(nny(j,k),3),num2str(SY(j,k)),'Color','k','FontWeight','bold','FontSize',20);
                end
            end
            
            
            axis equal;axis tight;title('X-Y Levenstein');
            drawnow
        end
    end
end


end