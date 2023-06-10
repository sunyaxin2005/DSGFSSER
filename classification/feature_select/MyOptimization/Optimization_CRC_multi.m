function [LTests] = Optimization_CRC_multi(XTrain, LTrain, XTest)
%OPTIMIZATION_CRC_MULTI 此处显示有关此函数的摘要
%   此处显示详细说明
UL = unique(LTrain);
K = length(UL);
NTset = size(XTest, 1);
LTest = zeros(NTset, 1);
[n, m] = size(XTrain);
X = XTrain';
Y = XTest';
XtX = X'* X;
XtY = X'* Y;
I = eye(n, n);
lamda = 0.5;
beta = 0.1;
T = XtX + lamda * I;
A = T^-1 * XtY;

num = 3;
LTests = cell(num, 1);
As = cell(num, 1);
A1 = A;
As{1} = A;
for jj = 2 : num
    for ii = 1 : size(Y, 2)
        Q  = 0;
        for kk = 1 : jj-1
            A = As{kk};
            Q = Q + A(:, ii) * A(:, ii)';
        end
        A1(:, ii) = (T + beta * Q)^-1* X'*Y(:, ii);
%         for kk = 1 : jj-1
%             A = As{kk};
%             Q = Q + A(:, ii);
%         end
%         A1(:, ii) = T* (X'*Y(:, ii) - Q);
    end
    As{jj} = A1;
end
for jj = 1 : num
    coefs = As{jj};
    coefsas1 = As{1};
    res = zeros(K,1);
    res1 = zeros(K,1);
    for ii = 1 : size(Y, 2)
        for kk = 1 : K
            coef1 = coefs(:, ii);
            coef1(LTrain ~= UL(kk)) = 0;
            
%             coef2 = coefsas1(:, ii);
%             coef2(LTrain ~= UL(kk)) = 0;           
            res(kk) = norm(Y(:, ii) - X * coef1,2);
            res1(kk) = norm( X * coef1,2);
        end
        [C,I] = min(res);
        LTest(ii) = UL(I);        
    end
    LTests{jj} = LTest;
end

end

