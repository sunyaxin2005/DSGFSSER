function [ UBMs] = GMM_UBM_TRAIN( features, train_inxs)
%GMM_UBM_TRAIN Summary of this function goes here
%   Detailed explanation goes here
%��Ҫѵ��GMM-UBM-MAP��������
X = features;
UBMs = cell(length(train_inxs), 1);
fprintf('GMM_UBM_TRAIN...');
for inx = 1 : length(train_inxs)
    fprintf('%d', inx);
    train = train_inxs{inx};
    Xtrain=X(train); 
    
    %��Ϊ�ǰ�������ѵ�����ݺͲ������ݵģ�������Ҫ�������ӳ�һ������
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
