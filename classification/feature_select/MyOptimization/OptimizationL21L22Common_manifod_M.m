function [Ws] = OptimizationL21L22Common_manifod_M( X, Y, L ,Para )
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
miu = 0.1;
lamda = Para.lambda;
beta = Para.beta;
%min_W,E ||X'W-Y||_2,1 + ||W||_2,1
%X:特征选择时， 样本个数*样本维数
%X稀疏表示时，样本维数*样本个数
%Y:
%X = X';
n = size(X, 1);
m = (size(X, 2));
I1 = gpuArray(eye(m, m));
Ut = gpuArray(eye(m, m));
Ut1 = gpuArray(eye(m, m));
XtY = gpuArray(X'* Y);
I = gpuArray(eye(n, n));
XtX = gpuArray(X'* (I + miu * L) * X) ;

Wi = zeros(size(XtY));
Ws = cell(3, 1);
for ii = 1 : 9
    for t=1:15
%         if(sum(sum(Wi)) ~= 0 )
%             W = (Ut * XtX+ 0.1 * Ut1)^-1 * (Ut * XtY);
%         else
%             W = (Ut * XtX+ lamda * I1 )^-1 * (Ut * XtY);
%         end
%         Ut = diag(sqrt(sum(W.^2, 2)).*2);
        
        W = (Ut * Ut1* XtX + lamda * Ut1 * I1)^-1 * (Ut * Ut1 * XtY - beta  * Ut * Wi);
        Ut = diag(sqrt(sum(W.^2, 2)).*2);
        if(sum(sum(Wi)) ~= 0 )
             Ut1 = diag(sqrt(sum(abs(W) .* abs(Wi), 2)).*2);
        end
%     Ut = Dt * A' * (A * Dt * A')^-1 * Y;
%     Dt = diag(sqrt(sum(Ut.^2, 2)).*2);
    end
    Wi = Wi + W; 
    %Wi = Wi./max(max(Wi));
    Ws{ii} = W;
    %WitWi = Wi*Wi';
    Ut = gpuArray(eye(m, m));
    Ut1 = gpuArray(eye(m, m));
    %aaa = sqrt(sum(Wi.^2, 2)).*2;
    %Ut1 = diag(aaa)./(sum(aaa));
end
end

