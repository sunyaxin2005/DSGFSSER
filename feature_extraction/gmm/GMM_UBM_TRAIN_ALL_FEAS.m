function [ UBM ] = GMM_UBM_TRAIN_ALL_FEAS( features )
%GMM_UBM_TRAIN_ALL_FEAS Summary of this function goes here
%   Detailed explanation goes here
X = features;
fprintf('GMM_UBM_TRAIN...');
for inx = 1 : 1
    fprintf('%d', inx);
    Xtrain=X; 
    
    %因为是按照语句分训练数据和测试数据的，所以需要将其连接成一个矩阵
    file_num = length(Xtrain);
    matrix_features = [];
    for file_inx = 1 : file_num
        feature = Xtrain{file_inx};
        feature(isnan(feature)) = 0;
        feature(isinf(feature)) = 0;
        matrix_features = [matrix_features;feature];
    end
    
    [ UBM] = GMM_UBM( matrix_features);
end

end

