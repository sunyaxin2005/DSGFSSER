function [ M_speakers ] = GMM_UBM_MAP( features, UBM )
%GMM_UBM_MAP Summary of this function goes here
%   Detailed explanation goes here
%%%%%%
% adapt UBM to suspect model   
avoid_non_psd = 1e-12;% avoid non-positive-semi-definite covariance matrices
max_adaption_iter = 15;    
file_num = length(features);
M_speakers = [];
for file_inx = 1 : file_num
    fprintf('%d ', file_inx)
    feature = features{file_inx};
    %feature = (feature - min(min(feature)))/(max(max(feature)) - min(min(feature)));
    feature(isnan(feature)) = 0;
    feature(isinf(feature)) = 0;
    M_speaker = UBM;    
    for Iadapt = 1:max_adaption_iter
        M_speaker = adaptUBM(M_speaker, feature, avoid_non_psd);
    end
    A.data = M_speaker;
    M_speakers = [M_speakers;A];
end

end

