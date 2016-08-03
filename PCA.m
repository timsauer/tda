function [X,s,U,mu] = PCA(Y,p)
%%% Algorithm 1.1 Principal Component Analysis (PCA)
%%% Inputs: Y - n-by-N matrix with columns y_i in R^n for i=1,...,N
%           p - percentage of variance to maintain (or dimension if p>1)
%%% Output: X - m-by-N matrix with columns x_i in R^m for i=1,...,N

mu = mean(Y,2);
Y = Y - repmat(mu,1,size(Y,2));  % Center the data
[U,S,~] = svd(Y,'econ');                       % Compute the SVD
S=diag(S).^2;                           % Extract variance of principal components
s=sum(S);                               % Total variance
s=cumsum(S/s);                          % Percent variance vs. number of components
if (p<1)
    m=find(s>=p,1);
else
    m=p;
end
X = U(:,1:m)'*Y;                        % Project onto first m principal components
