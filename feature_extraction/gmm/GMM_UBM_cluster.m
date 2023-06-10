function [ feature_cluster ] = GMM_UBM_cluster( features, UBM )
%GMM_UBM_CLUSTER Summary of this function goes here
%   Detailed explanation goes here
weights = UBM.PComponents;
[D, INX] = sort(weights, 'descend');
mu =  UBM.mu;
cluster_center = mu(INX(1:10), :);

files_num = length(features);
feature_cluster = cell(files_num, 10);
for file_inx = 1 : files_num
    feature = features{file_inx};
    feature = (feature - min(min(feature)))/(max(max(feature)) - min(min(feature))); 
    [frame_num, feature_dim] = size(feature);
    pds = pdist2(feature, cluster_center, 'cosine');
    pds = sort(pds');
    pds = pds';
end

end

