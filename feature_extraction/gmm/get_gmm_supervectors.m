function [ stats_features,  M_speakers] = get_gmm_supervectors(features, UBMS)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%传进来的数据类型为CELL，每一个CELL内的数据为：数据的个数*数据的维数
M_speakers = cell(length(UBMS), 1);
stats_features = cell(length(UBMS), 1);
for inx = 1 : length(UBMS)
    UBM = UBMS{inx};
    M_speaker = GMM_UBM_MAP( features, UBM );
    stats_feature = gmm_supervector( M_speaker); 
    
    stats_feature(isnan(stats_feature)) = 0;
    stats_feature(isinf(stats_feature)) = 0;
    M_speakers{inx} = M_speaker;
    stats_features{inx} = stats_feature;   
end
end

