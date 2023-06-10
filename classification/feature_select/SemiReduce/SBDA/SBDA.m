function [W,C]=SBDA(Sb,Sw,St,options)
%sparse undirect learning framework/undirect sparse projection learning
%framework:SULF/USPL

%min miu*w^T*Sw*w+(1-miu)*(c-w)^T*Sw*(c-w)+alfa*||c||^2+beta*|c|
%s.t.      w^T*Sb*w=1

%X:column vectors as a pattern to form the matrix X
%card: the number of non-zero elements in the projection C or W
%OUTPUT: C ,which is used for dimensionality reduction as PCA projection
%Please see the reference:Zhihui Lai, Yong Xu, IEEE, Zhong Jin, and David Zhang,Human Gait Recognition via Sparse
%Discriminant Projection Learning,IEEE transactions on Circuits and Systems
%for Video Technology,2014

miu=options.miu;card=options.card;lambda=options.lambda;LagMutiplier=options.LagMutiplier;ite=options.ite;
if miu>=1
    error('wrong parameters')
end

% if Sb==Sw || Sw==St
%     Sb=M;Sw=M;St=eye(size(Sb));
% end
[U D V]=svd(Sw);%select a better initialization value for C and save time in iteration
C=U;%used for initialization C
D=diag(sqrt(diag(D)));%preparisions for iterations
flag=0;rss=[];%LagMutiplier=0.001;%initailization
alfa=1000;Wold=[];
while 1
    %fix C to get W
    W=(miu*Sb+(1-miu)*Sw+LagMutiplier*St+0.00001*eye(size(Sw)))\((1-miu)*Sw*C);
    %fix W to get C
    Y=(U*D)'*W;
%     Y=Y./(ones(size(Y,1),1)*sqrt(sum(Y.*Y)));
    X=(U*D)';
%     X=X./(ones(size(X,1),1)*sqrt(sum(X.*X)));
    for i=1:size(U*D)
        B = larsen(X,Y(:,i),lambda,-card);%inf
        C(:,i) = B(end,:)'/norm(B(end,:));
    end 
    RSS = norm(Y-X*C,2)^2 + lambda*norm(C,2)^2+sum(sum(abs(C)));
    rss(ite) = RSS;
    normW(ite)=norm(W,2);
    
    if ite>1
        if abs((rss(ite)-rss(ite-1))/rss(ite))<0.001%|| norm(W-Wold,2)<0.01
            flag=1;
            break
        end
    end
    if ite>100
        break
    end
    ite=ite+1;
    Wold=W;
end
% if flag==1
%     disp('auto convergence and break')
% else
%     disp('achieve the number of iterations setting: 100. forced break')
% end
% 1
    
    