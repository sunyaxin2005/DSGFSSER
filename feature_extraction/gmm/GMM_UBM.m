function [ UBM] = GMM_UBM( matrix_features)
%GMM_UBM Summary of this function goes here
%   Detailed explanation goes here
%%%%%%
%��һ��������UBM
%GMM_UBM��һЩԤ�����
num_gaussians = 32;
avoid_non_psd = 1e-12;% avoid non-positive-semi-definite covariance matrices
stat_options = statset('MaxIter', 500);
max_adaption_iter = 15;%��MAP���ʱִ��15��

%����UBM
lastwarn('');
UBM = gmdistribution.fit(matrix_features, num_gaussians, 'CovType', 'diagonal', 'Regularize', avoid_non_psd, 'Options', stat_options);
if ~isempty(lastwarn), beep, return, end % STOP here is there is a warning (failure to converge)

end