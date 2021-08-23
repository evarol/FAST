function [bhat,P,sbest,q]=rrwoc_wahba(x,y,delta,nu,k,max_iter)
% x     - source
% y     - target
% delta - probability of success
% nu    - margin
% k     - expected number of outliers (for assessing randomized iterations)
warning('off');
[n,d]=size(y);
[m,~]=size(x);

iter=ceil(log(1-delta)/log(1-(prod(m-k:-1:m-k-d+1)/(prod(m:-1:m-d+1)*prod(n:-1:n-d+1)))));
sbest = 0;
inlier = [];
for q = 1: min(max_iter, iter)
    if mod(q, 10) == 0
        disp([num2str(q) '/' num2str(iter) '-- best score: ' num2str(sbest)]);
    end
    
    %% candidate selection
    
%     CX = randperm(size(x,1)); CX = CX(1:d+1);
r=randi(size(x,1));
[~,idx]=sort(pdist2(x(r,:),x),'ascend');
CX=idx(1:d+1);
    %     CY = randperm(size(y,1)); CY = CY(1:d+1);
    %
    %     [R,T] = wahba(x(CX,:),y(CY,:));
    xsubset=x(CX,:);
    ysubset=volumetric_sampling(xsubset,y);
    
    %% Iterate to find a best subset possible
    
    for t=1:4
        for i=1:size(ysubset,1)
            [ysubset_candidate{i},hardGW_candidate(i,1)]=volumetric_sampling(xsubset,setdiff(y,ysubset(i,:),'rows'),ysubset(i,:));
        end
        [hardGW,idx]=min(hardGW_candidate);
        ysubset=ysubset_candidate{idx};
    end
     
    
    [R,T] = wahba(xsubset,ysubset);
    
    yhat = x*R + T;
    D = pdist2(y,yhat);
    s = sum(any(min(D, [], 2) <= nu, 2));
    if s > sbest
%         subplot(1,2,1)
%         cla
%         hold on
%         plot3(yhat(CX,1),yhat(CX,2),yhat(CX,3),'b.','MarkerSize',20);
%         plot3(ysubset(:,1),ysubset(:,2),ysubset(:,3),'r.','MarkerSize',20);
%         legend('moved','target')
%         set(gca,'FontWeight','bold','FontSize',20,'TickLength',[0 0]);set(gcf,'Color','w');axis equal;axis tight
%         subplot(1,2,2)
%         cla
%         hold on
%         plot3(yhat(:,1),yhat(:,2),yhat(:,3),'b.','MarkerSize',20);
%         plot3(y(:,1),y(:,2),y(:,3),'r.','MarkerSize',20);
%         legend('moved','target')
%         set(gca,'FontWeight','bold','FontSize',20,'TickLength',[0 0]);set(gcf,'Color','w');axis equal;axis tight
%         drawnow
        [P.matched_cells_y,P.matched_cells_x]=find(D<= nu);
        sbest = s;
        bhat=[R;T];
    end
    if sbest >= m-k
        break;
    end
end
disp(['Algorithm reached an inlier set after ' num2str(q) ' iterations']);

% [~,~,bhat]=wahba(Pinlier(inlier,:)*y(:,1:end-1),x(inlier,:));
% inlier_set = inlier;
inlier_set=[];


end
