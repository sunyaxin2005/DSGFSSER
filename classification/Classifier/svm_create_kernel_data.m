function [ Ktest ] = svm_create_kernel_data( Xtrain, Xtest, kernel_type)
%SVM_CREATE_KERNEL_DATA Summary of this function goes here
%   Detailed explanation goes here
%Xtrain:sample_num * feature_dim
%Xtest : sample_num * feature_dim
test_num = size(Xtest, 1);
train_num = size(Xtrain, 1);
switch kernel_type
    case 'linear'
         Ktest = Xtest * Xtrain';
    case 'example1'
          % 使用的核函数 K(x,x') = ||x|| * ||x'||
        Ktest = ones(test_num,train_num);
        for i = 1:test_num
            for j = 1:train_num
                Ktest(i,j) = sum(Xtest(i,:).^2)^0.5 * sum(Xtrain(j,:).^2)^0.5;
            end
        end
        Ktest = [(1:size(Ktest, 1))', Ktest];
    case 'example2'
         % 使用的核函数 K(x,x') = (x * x') / ||x|| * ||x'|| 
        Ktest = ones(test_num,train_num);
        for i = 1:test_num
            for j = 1:train_num
                Ktest(i,j) = Xtest(i,:)*Xtrain(j,:)'/(sum(Xtest(i,:).^2)^0.5 * sum(Xtrain(j,:).^2)^0.5);
            end
        end
        Ktest = [(1:size(Ktest, 1))', Ktest];
    case 'poly_nomial'
        Ktest = (Xtest*Xtrain' + 1).^2;
        Ktest = [(1:size(Ktest, 1))', Ktest];
    case 'KL'
         Ktest = ones(test_num,train_num);
        %映射训练数据的核函数
        for i = 1:test_num
            wi = Xtest(i).data.PComponents';
            wi = repmat(wi, 1, size(Xtrain(i).data.mu, 2));
            mi = Xtest(i).data.mu;
            signai = zeros(size(wi));
            signai = signai';
            signai(:, :) =  Xtest(i).data.Sigma(1, :, :);
            signai = signai';    
            vi = mi./ sqrt(signai) .* wi;
            for j = 1:train_num
                wj = Xtrain(j).data.PComponents';
                wj = repmat(wj, 1, size(Xtrain(j).data.mu, 2));
                mj = Xtrain(j).data.mu;
                signaj = zeros(size(mj));
                signaj = signaj';       
                signaj(:, :) =  Xtrain(j).data.Sigma(1, :, :);
                signaj = signaj';  
                vj = mj./ sqrt(signaj) .* wj;
                Ktest(i,j) = sum(sum(vi.*vj));
            end
        end
        Ktest = [(1:size(ktest2, 1))', ktest2];
    case 'HIK'
        %Ktest = ones(test_num,train_num);
        Ktest = zeros(test_num, train_num);
        for i = 1:test_num
            for j = 1:train_num
                for inx = 1 : size(Xtrain, 2)
                    Ktest(i, j) = Ktest(i, j) + min(Xtest(i,inx).^2, Xtrain(j,inx).^2);
                end
            end
        end
        %Ktest = [(1:size(Ktest, 1))', Ktest];
end

end

