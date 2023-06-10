function [ UBMs] = GMM_UBM_TRAIN( features, train_inxs)
%GMM_UBM_TRAIN Summary of this function goes here
%   Detailed explanation goes here
%需要训练GMM-UBM-MAP的特征，
X = features;
UBMs = cell(length(train_inxs), 1);
fprintf('GMM_UBM_TRAIN...');
for inx = 1 : length(train_inxs)
    fprintf('%d', inx);
    train = train_inxs{inx};
    Xtrain=X(train); 
    
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
    UBMs{inx} = UBM;
end
end

