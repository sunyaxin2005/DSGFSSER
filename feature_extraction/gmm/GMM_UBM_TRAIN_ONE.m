function [ UBM ] = GMM_UBM_TRAIN_ONE( Xtrain, gnd, ii,files_people_label )
%GMM_UBM_TRAIN_ONE Summary of this function goes here
%   Detailed explanation goes here
fprintf('GMM_UBM_TRAIN_ONE.....');
file_num = length(Xtrain);
feature = Xtrain{1};
matrix_features = zeros(1000000, size(feature, 2));
UL = unique(gnd);
start_I = 1;
for file_inx = 1 : file_num
    if gnd(file_inx) == UL(ii) && files_people_label(file_inx) > 25
        feature = Xtrain{file_inx};
        % feature = (feature - min(min(feature)))/(max(max(feature)) - min(min(feature)));  
        feature(isnan(feature)) = 0;
        feature(isinf(feature)) = 0;
        end_I = start_I + size(feature, 1) - 1;
        matrix_features(start_I:end_I, :) = feature;
        start_I = end_I + 1;
    end
end
matrix_features = matrix_features(1:end_I, :);
[ UBM] = GMM_UBM( matrix_features);
end

