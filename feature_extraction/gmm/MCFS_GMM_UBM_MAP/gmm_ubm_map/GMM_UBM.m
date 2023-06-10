function [ UBM] = GMM_UBM( matrix_features)
%GMM_UBM Summary of this function goes here
%   Detailed explanation goes here
%%%%%%
%第一步：建立UBM
%GMM_UBM的一些预设参数
num_gaussians = 32;
avoid_non_psd = 1e-12;% avoid non-positive-semi-definite covariance matrices
stat_options = statset('MaxIter', 500);
max_adaption_iter = 15;%在MAP拟合时执行15次

%生成UBM
lastwarn('');
UBM = gmdistribution.fit(matrix_features, num_gaussians, 'CovType', 'diagonal', 'Regularize', avoid_non_psd, 'Options', stat_options);
if ~isempty(lastwarn), beep, return, end % STOP here is there is a warning (failure to converge)

end