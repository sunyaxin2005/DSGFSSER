function [ W] = getrelation_eig(  XL, L, XUL )
%LOCATELRGPSSDR Summary of this function goes here
%   Detailed explanation goes here
[eigvector, W] = semi_dimension_reduction_main( XL, L, XUL, 'LDA');
fea1 = [XL;XUL];
W = fea1 * eigvector * W;
% XL = XL * eigvector * W(:, 1 : 100);
% XUL = XUL * eigvector * W(:, 1 : 100);
% fea1 = [XL;XUL];
% distance = pdist2(fea1, fea1);
% [sorted,index] = sort(distance);
% sita = sorted(5+1, :);
% sita = sita' *  sita;
% W = exp(-(distance.^2)./sita);

end

