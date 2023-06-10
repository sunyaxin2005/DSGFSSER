function [ X ] = sun_norm_by_l2norm( X )
%SUN_NORM1 Summary of this function goes here
%   Detailed explanation goes here
%X: N¡ÁD
for ii = 1 : size(X, 1)
    X(ii, :) = X(ii, :)./norm(X(ii, :), 2);
end
% maxd = max(X);
% mind = min(X);
% maxd = repmat(maxd, size(X, 1), 1);
% mind = repmat(mind, size(X, 1), 1);
% X = (X - maxd)./(maxd - mind);
% X(isnan(X)) = 0;
% X(isinf(X)) = 0;

end

