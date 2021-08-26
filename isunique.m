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