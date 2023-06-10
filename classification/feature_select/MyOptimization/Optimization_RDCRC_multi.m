function [LTests] = Optimization_RDCRC_multi(XTrain, LTrain, XTest)
%OPTIMIZATION_CRC_MULTI 此处显示有关此函数的摘要
%   此处显示详细说明
gama = 0.5;
XTrain = XTrain';
train_num = size(XTrain, 2);
ex_data1 = [];
UL = unique(LTrain);
class_num = length(UL);
M=eye(train_num);
start_inx = 1;
L_new = LTrain;
for ii=1:class_num
    xi = XTrain(:, LTrain == UL(ii));
    Aaa = xi'*xi;
    dim = size(Aaa, 1);
    M(start_inx:start_inx + dim - 1,start_inx:start_inx + dim - 1) = Aaa;
    L_new(start_inx:start_inx + dim - 1) = UL(ii);
    start_inx = start_inx + dim;
    ex_data1 = [ex_data1, xi];
end
XTrain = ex_data1';
LTrain = L_new;

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
lamda = 0.05;

T = XtX + lamda * I;% + gama * XtX- gama*M * 2;
A = T^-1 * XtY;
beta = trace(T);
num = 3;

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
        Q = Q/trace(Q);
        A1(:, ii) = (T + beta * Q)^-1* X'*Y(:, ii);
%         for kk = 1 : jj-1
%             A = As{kk};
%             Q = Q + A(:, ii);
%         end
%         A1(:, ii) = T* (X'*Y(:, ii) - Q);
    end
    As{jj} = A1;
end

LErr = zeros(NTset, 1);
LErrs = zeros(NTset, num);
LTests = zeros(NTset, num);
res_ALL = zeros(size(Y, 2), K, num);
for jj = 1 : num
    coefs = As{jj};
    %coefsas1 = As{1};
    res = zeros(K,1);
    res1 = zeros(K,1);
    for ii = 1 : size(Y, 2)
        for kk = 1 : K
            coef1 = coefs(:, ii);
            coef1(LTrain ~= UL(kk)) = 0;
            
%             coef2 = coefsas1(:, ii);
%             coef2(LTrain ~= UL(kk)) = 0;           
            res(kk) = norm(Y(:, ii) - X * coef1,2);
            res_ALL(ii, kk, jj) = res(kk);
            res1(kk) = norm( X * coef1,2);
        end
        [C,I] = min(res);
        LTest(ii) = UL(I);   
        LErr(ii) = C;
    end
    LTests(:, jj) = LTest;
    LErrs(:, jj) = LErr;
end
res_ALL1 = zeros(size(Y, 2), K);
res_ALL1 = sum(res_ALL,3); 
end
