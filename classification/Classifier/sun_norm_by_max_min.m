function [ X ] = sun_norm_by_max_min( X )
%SUN_NORM1 Summary of this function goes here
%   Detailed explanation goes here
%X: N¡ÁD
X(isnan(X)) = 0;
maxd = max(X);
mind = min(X);
inx = ((maxd -mind)  ~= 0 & ~isnan(maxd) & ~isinf(maxd) & ~isnan(mind) & ~isinf(mind));
maxd = maxd(inx);
mind = mind(inx);
X = X(:, inx);
maxd = repmat(maxd, size(X, 1), 1);
mind = repmat(mind, size(X, 1), 1);
X = (X - maxd)./(maxd - mind);
X(isnan(X)) = 0;
X(isinf(X)) = 0;

end

