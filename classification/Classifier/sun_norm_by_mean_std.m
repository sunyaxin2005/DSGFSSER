function [ X ] = sun_norm_by_mean_std( X )
%SUN_NORM Summary of this function goes here
%   Detailed explanation goes here
%传进来的数据为:N*D
X(isnan(X)) = 0;
stdd = std(X);
inx = (stdd ~= 0 & ~isnan(stdd) & ~isinf(stdd));
stdd = stdd(inx);
meandd = mean(X);
meandd = meandd(inx);
X = X(:, inx);
stdd = repmat(stdd, size(X, 1), 1);
meandd = repmat(meandd, size(X, 1), 1);
X = (X - meandd)./stdd;

end

