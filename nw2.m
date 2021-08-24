function [T,P,skip,xhat,match]=nw2(x,y,indelcost)
%needleman-wunsch algorithm to calculate levenstein distance
% match cost = -|x-y|
% indel cost = -1


%% encode
T=nan(size(x,1)+1,size(y,1)+1);
P=nan(size(x,1)+1,size(y,1)+1,2);
skip=nan(size(x,1)+1,size(y,1)+1);
T(:,1)=indelcost*(0:-1:-size(x,1))';
T(1,:)=indelcost*(0:-1:-size(y,1));
for i=2:size(T,1)
    for j=2:size(T,2)
        [T(i,j),idx]=max([T(i,j-1)-indelcost,T(i-1,j)-indelcost,T(i-1,j-1)-abs(x(i-1)-y(j-1))]); %%cost could be how parallel we are
        if idx==1
            P(i,j,1)=i;
            P(i,j,2)=j-1;
            skip(i,j)=1;
        elseif idx==2
            P(i,j,1)=i-1;
            P(i,j,2)=j;
            skip(i,j)=1;
        elseif idx==3
            P(i,j,1)=i-1;
            P(i,j,2)=j-1;
            skip(i,j)=0;
        end
    end
end

%% decode
% currentpoint=[size(T,1) size(T,2)];
[xmin,xminpos]=max(T(end,:));
[ymin,yminpos]=max(T(:,end));
if xmin>ymin
    currentpoint=[size(T,1) xminpos];
else
    currentpoint=[yminpos size(T,2)];
end
xhat=nan(size(x));
match=[];
while and(currentpoint(1)~=1,currentpoint(2)~=1)
    if skip(currentpoint(1),currentpoint(2))~=1
        xhat(currentpoint(1)-1)=y(currentpoint(2)-1);
        match=[match;currentpoint(1)-1 currentpoint(2)-1];
    end
    currentpoint=squeeze(P(currentpoint(1),currentpoint(2),:));
    if or(currentpoint(1)==1,currentpoint(2)==1)
        break
    end
end

return
% 
% beta=nanmedian(x-xhat);
% 
% close all
% hold on
% imagesc(T)
% currentpoint=[size(T,1)-1 size(T,2)-1];
% while and(currentpoint(1)~=1,currentpoint(2)~=1)
%     plot(currentpoint(1),currentpoint(2),'r.')
%     currentpoint=squeeze(P(currentpoint(1),currentpoint(2),:));
%     drawnow
%     if or(currentpoint(1)==1,currentpoint(2)==1)
%         break
%     end
% end
end
