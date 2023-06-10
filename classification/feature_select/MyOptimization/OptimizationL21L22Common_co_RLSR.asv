function [ W2 ] = OptimizationL21L22Common_co_RLSR( XL, YL1, YL2, XU, YU1, YU2, p,gamma, MAX_ITER)
%OPTIMIZATIONL21L22COMMON Summary of this function goes here
%   Detailed explanation goes here

epsilon=1e-5;
if nargin<7
    p=1;
end

if nargin<8
    gamma=1;
end

if nargin<9
    MAX_ITER=100;
end

q = 2/p - 1;

nl=size(XL,2);
num=nl+size(XU,2);
d=size(XL,1);
c=size(YL,2);
X=zeros(size(XL,1),num);
X(:,1:nl)=XL;
X(:,nl+1:num)=XU;

Y1=ones(num,c)/c;
Y1(1:nl,:)=YL1;
Y1(nl+1:num,:)=YU1;

Y2=ones(num,c)/c;
Y2(1:nl,:)=YL2;
Y2(nl+1:num,:)=YU2;

H = eye(num) - ones(num)/num;
bigTheta=eye(d)/d;

obj=zeros(MAX_ITER,1);

XHX=X*H*X';
% 记住了XHX和XHY是变量名。
for iter=1:MAX_ITER
    XHY1=X*H*Y1;
    W1=(XHX+gamma*(bigTheta^-2))\XHY1;
    temp=sum(W1.*W1,2).^(1/(q+1))+epsilon;
    bigTheta=diag( temp/ sum(temp) ).^(q/2);
    b=(sum(Y1,1)'-sum(W1'*X,2))/num;

    % updata Yu
    for i=nl+1:num
        Y1(i,:)=X(:,i)'*W1+b';
        Y1(i,:) = EProjSimplex_new(Y1(i,:));
    end
    obj(iter)=F22norm(X'*W1+repmat(b',[num 1])-Y1)+gamma*L2Pnorm(W1,p)^2;
    
    if iter ==1
        minObj=obj(iter);
        bestW1 = W1;
    else
        if ~isnan(obj(iter)) && obj(iter) <= minObj
            minObj=obj(iter);
            bestW1 = W1;
        end
    end
    
    if iter>1
        change=abs((obj(iter)-obj(iter-1))/obj(iter));
        if change<1e-8
            break;
        end
    end
end

W=bestW;

theta=sum(W.*W,2).^(p/2);
theta=theta/sum(theta);
[~, ranked] = sort(theta, 'descend');
end

