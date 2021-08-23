function [vol_dist,P]=volumetric_distance_2(X,Y,numperm)
if nargin<3
    numperm=2;
end
[~,xnn]=sort(squareform(pdist(X)),2,'ascend');
vol_dist=zeros(size(X,1),size(Y,1));
bestsofar=Inf;
P=nan(size(X,1),size(Y,1),numperm);
for t=1:numperm
    for i=1:size(X,1)
        Xsubset=X(xnn(i,[1 randsample(2:10,3)]),:);
        for j=1:size(Y,1)
            [Ysubset,hardGW]=volumetric_sampling(Xsubset,Y,j);
            for k=1:size(Ysubset,1)
                idxy=find(norms(Y-Ysubset(k,:),2)==0);
                idxx=find(norms(X-Xsubset(k,:),2)==0);
                P(idxx,idxy,t)=norm(Ysubset(k,:)-Xsubset(k,:));
            end
            [R,T,beta]=wahba(Xsubset,Ysubset);
            Xsubsethat=Xsubset*R + T;
            vol_dist(i,j,2,t)=norm(Ysubset-Xsubsethat,'fro');
            if vol_dist(i,j,2)<bestsofar;
                bestsofar=vol_dist(i,j,2,t);
                subplot(1,2,2);
                cla
                hold on
                plot3(Xsubsethat(:,1),Xsubsethat(:,2),Xsubsethat(:,3),'b.','MarkerSize',20);
                plot3(Ysubset(:,1),Ysubset(:,2),Ysubset(:,3),'r.','MarkerSize',20);
                set(gca,'Color','k');
                axis equal;axis tight;axis off;title([num2str([i j]) ' average distance: ' num2str(bestsofar)]);
                subplot(1,2,1);
                Xhat=X*R+T;
                cla
                hold on
                plot3(Xhat(:,1),Xhat(:,2),Xhat(:,3),'b.','MarkerSize',20);
                plot3(Y(:,1),Y(:,2),Y(:,3),'r.','MarkerSize',20);
                set(gca,'Color','k');
                axis equal;axis tight;axis off;title([num2str([i j]) ' average distance: ' num2str(bestsofar)]);drawnow
            end
            vol_dist(i,j,1,t)=hardGW;
            [i j]
        end
    end
end
end

function n=norms(x,dim)
n=sqrt(sum(x.^2,dim));
end
