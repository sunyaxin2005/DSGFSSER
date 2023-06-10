function [ classes, accuracy, dec_values] = svm_using_kernel_test( Xtest, Ltest,Xtrain, model, kernel_type )
%SVM_USING_KERNEL_TEST Summary of this function goes here
%   Detailed explanation goes here
%Xtrain:sample_num * feature_dim
%Xtest : sample_num * feature_dim
[ Ktest ] = svm_create_kernel_data( Xtrain, Xtest, kernel_type);
[classes, accuracy, dec_values] = svmpredict(Ltest, Ktest, model);  
end

