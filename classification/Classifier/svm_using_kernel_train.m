function [ model_precomputed ] = svm_using_kernel_train( Xtrain, Ltrain, kernel_type )
%SVM_USING_KERNEL_TRAIN Summary of this function goes here
%   Detailed explanation goes here
%Xtrain:sample_num * feature_dim
%Xtest : sample_num * feature_dim
[ Ktrain ] = svm_create_kernel_data( Xtrain, Xtrain, kernel_type);
model_precomputed = svmtrain(Ltrain, Ktrain, '-c 16 -t 4');
end

