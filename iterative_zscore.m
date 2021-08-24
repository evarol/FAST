function X=iterative_zscore(X)

for t=1:100
    X=zscore(X,[],1);
    X=zscore(X,[],2);
end