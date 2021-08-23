function [R,T,beta]=weighted_wahba(X,Y,w)
% returns solution to ||Y - X*R - T||, R \in S0(3)
% alternate form ||Y - [X 1]*beta||


X0=X-mean(X,1);
Y0=Y-mean(Y,1);
B=0;
for i=1:size(X,1)
    B=w(i)*X0(i,:)'*Y0(i,:);
end
[U,~,V]=svd(B);

M=[1 0 0;0 1 0;0 0 det(U)*det(V)];
R=U*M*V';


T= mean(Y - X*R,1);

beta = [R;T];
end