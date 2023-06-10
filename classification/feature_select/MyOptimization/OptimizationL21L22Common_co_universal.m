function [ W, Ws ] = OptimizationL21L22Common_co_universal( X, Ys, Para)
%OPTIMIZATIONL21L22COMMON Summary of this function goes here
%   Detailed explanation goes here
miu = 0.1;
lamda = Para.lambda;
beta = Para.beta;
%min_W,E ||X'W-Y||_2,1 + ||W||_2,1
%X:特征选择时， 样本个数*样本维数
%X稀疏表示时，样本维数*样本个数
%Y:
%X = X';
n = size(X, 1);
m = size(X, 2);
I = eye(n, n);
Uts = cell(length(Ys), 1);
for ii = 1 : length(Ys)
    Uts{ii} = eye(m, m);
end
XtYs = cell(length(Ys), 1);
for ii = 1 : length(Ys)
   XtY =X'* Ys{ii};
   XtYs{ii} = XtY;
end
XtX = X'* X ;
Ut = eye(m, m);
I1 = eye(m, m);
Ws =  cell(length(Ys), 1);
for t=1:30
    for ii = 1 : length(Ys)
        Ws{ii} =  (Uts{ii}*Ut * XtX+ lamda * Uts{ii} + beta * Ut + 0.00000001 * I1)^-1 * Uts{ii}*Ut *XtYs{ii};
        Uts{ii} = diag(sqrt(sum(Ws{ii}.^2, 2)).*2);   
    end
    W = [];
    for ii = 1 : length(Ys)
        W = [W, Ws{ii}];
    end
    Ut = diag(sqrt(sum(W.^2, 2)).*2);
end
end

