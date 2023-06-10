function [ classes ] = svm_test_main( Xtest, Ltest, Xtrain, model,  kernel_type)
%SVM_TEST_MAIN Summary of this function goes here
%   Detailed explanation goes here
switch kernel_type
    case 'rbf'
        [classes, accuracy, dec_values] = svmpredict(Ltest, Xtest, model);  
    case {'linear', 'example1', 'example2', 'poly_nomial', 'KL', 'HIK'}
        [ classes, accuracy, dec_values] = svm_using_kernel_test( Xtest, Ltest,Xtrain, model, kernel_type );
end
end

