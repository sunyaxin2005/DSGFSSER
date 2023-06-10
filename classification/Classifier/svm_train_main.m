function [ model,  train_param] = svm_train_main( Xtrain, Ltrain, kernel_type )
%SVM_TRAIN_WITH_SEL_PARAM Summary of this function goes here
%   Detailed explanation goes here
switch kernel_type
    case 'rbf'
        clear global gTrainingFitnessIndex
         global gTrainingFitnessIndex;global gParms;
         for ntimes=1:1
             kkfold=5;
             tindices = crossvalind('Kfold',Ltrain,kkfold);
             for jj=1:kkfold
                 ttest = (tindices == jj); ttrain = ~ttest;       
                 gTrainingFitnessIndex(:,jj+(ntimes-1)*kkfold)=[find(ttrain)' find(ttest)']; %randperm(length(Ltrain));
             end;
         end
         methodid=66;
         train_params=select_params_bycv_for_lrt_lmc(Xtrain',Ltrain,methodid);
         cc=train_params(1);
         gamm=train_params(2);   
%         cc = 6;
%         gamm = -6;
         kerneltype=2;%cc=100;gamm=0.5;
         train_param = sprintf('-c %d -g %d -t %d', 2^cc,2^gamm,kerneltype);
        %[bestacc,bestc,bestg] = svm_cg(Ltrain,Xtrain);
        %train_param = sprintf('-c %d -g %d -t %d', bestc, bestg, 2);
        model = svmtrain(Ltrain,Xtrain,train_param);
    case {'linear', 'example1', 'example2', 'poly_nomial', 'KL', 'HIK'}
        train_param = '-t 4';
        [ model ] = svm_using_kernel_train( Xtrain, Ltrain, kernel_type );
end

end

