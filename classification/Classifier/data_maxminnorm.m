function [Xnorm]=data_maxminnorm(X)
%传进来的数据必须是：sample_num * dimension
X(isnan(X)) = 0;
X(isinf(X)) = 0;
tol=1e-3;
maxv = max(X);
minv = min(X);
maxv = repmat(maxv, size(X, 1), 1);
minv = repmat(minv, size(X, 1), 1);
Xnorm = (X - minv + tol)./(maxv-minv+tol);
end