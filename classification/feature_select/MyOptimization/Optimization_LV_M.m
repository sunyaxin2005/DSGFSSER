function [W] = Optimization_LV_M(X,Y, Para)
%OPTIMIZATIONNSSRD 此处显示有关此函数的摘要
%   此处显示详细说明

k=2;

%Ws = Get_01Weight( X, Para.gnd );
Ws = Get_01_manifold_weight( X, Para.gnd );
Ds = diag(sum(Ws(:,:),2));
Ls = (Ds)-Ws;

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
I1 = eye(m, m);
Ut1 = eye(m, m);
Ut2 = eye(m, m);
Ut3 = eye(m, m);
XtY = X'* Y;
XtX = X'* X ;
V = Y;
Vs = cell(k, 1);
for ii = 1 : k
    Vs{ii} = Y;
end

W1s = cell(k, 1);
Ut1s = cell(k, 1);

for ii = 1 : k
    Ut1s{ii} = Ut1;
end
WR =[];
for t=1:15
    for tt = 1 : k
        V = Vs{tt};
        VstVs = Y*Y';
        for tt1 = 1:k
            Vm = Vs{tt1};
            if tt1 ~= tt
                VstVs = VstVs + Vm*Vm';
            end
        end
        VstVs = VstVs./k;
        Ut1 = Ut1s{tt};
        W1 = (Ut1*Ut3 * XtX+ lamda * Ut1 + beta * Ut3 + 0.000001 * I1)^-1 * Ut1*Ut3 * X'*V;
        Ut1 = diag(sqrt(sum(W1.^2, 2)).*2);
        Ut1s{tt} = Ut1;
        W1s{tt} = W1;
        V = (VstVs+I+Ls)^-1*X*W1;
        V = V *sum(sum(abs(Y)))/sum(sum(abs(V)));
        Vs{tt} = V;
    end
%     if t== 1
%         Ut3 = Ut1;
%     else
%         Ut3 = diag(sqrt(sum(W.^2, 2)).*2);
%     end
    W2 = (Ut2*Ut3 * XtX+ lamda * Ut2 + beta * Ut3 + 0.000001 * I1)^-1 * Ut2*Ut3 * XtY;
    Ut2 = diag(sqrt(sum(W2.^2, 2)).*2);    
    %W = [W1, W2];    
    W = W2;
    for ii = 1 : k
        W1 = W1s{ii};
        W = [W1, W];
    end
    Ut3 = diag(sqrt(sum(W.^2, 2)).*2)/k;
%     W = (Ut* XtX+ lamda * I1)^-1 * Ut * XtY;
%     Ut = diag(sqrt(sum(W.^2, 2)).*2);
%     Ut = Dt * A' * (A * Dt * A')^-1 * Y;
%     Dt = diag(sqrt(sum(Ut.^2, 2)).*2);
    WR = W;
end
W = WR;
end

