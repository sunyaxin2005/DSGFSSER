function [ output_args ] = lgssdr( input_args )
%lgssdr Summary of this function goes here
%   Detailed explanation goes here

%generatesamples;

load yaleB.mat;

for kkk=0:5
    X=x';%训练集
    XT=xt';%测试集
    N=size(y,1);%训练数据个数
    K=kkk;%最近邻接点个数(不包括自身)
    dim=20;%目标维数
    NOC=20;%约束点对个数

    run=10;%实验次数
    myprecision=0;
    zhangprecision=0;

    options.PCARatio = 0.98;
    [eigvector, eigvalue] = PCA(X', options);
    X = X'*eigvector;
    XT=XT'*eigvector;
    X=X';
    XT=XT';

    S=eye(N,N);%邻接矩阵
    S=sparse(S);
    A=ones(N,N)-eye(N,N);%非邻接矩阵
    A=sparse(A);
    X2 = sum(X.^2,1);
    distance = repmat(X2,N,1)+repmat(X2',1,N)-2*X'*X;
    [sorted,index] = sort(distance);
    neighborhood = index(2:(1+K),:);
    for i=1:N
        for j=2:1+K
            S(i,index(j,i))=1;
            S(index(j,i),i)=1;
            A(i,index(j,i))=0;
            A(index(j,i),i)=0;
        end
    end

    for round=1:run
        S1=S;
        A1=A;
        %M=zeros(N,N);%must-link矩阵
        M=sparse(N,N);
        %C=zeros(N,N);%cannot-link矩阵
        C=sparse(N,N);
        nm=0;%must-link点对个数
        nc=0;%cannot-link点对个数

        for i=1:NOC
            part1=randint(1,1,[1,N]);
            part2=randint(1,1,[1,N]);
            while part1==part2 || M(part1,part2)==1 || C(part1,part2)==1
                part2=randint(1,1,[1,N]);
            end
            if y(part1)==y(part2)
                M(part1,part2)=1;
                M(part2,part1)=1;
                nm=nm+1;
            end
            if y(part1)~=y(part2)
                C(part1,part2)=1;
                C(part2,part1)=1;
                nc=nc+1;
            end
        end

        Ds = sum(S1(:,:),2);
        Ls = diag(Ds)-S1;
        Dm = sum(M(:,:),2);
        Lm = diag(Dm)-M;
        Dc = sum(C(:,:),2);
        Lc = diag(Dc)-C;
        Da = sum(A1(:,:),2);
        La = diag(Da)-A1;

        A2=X*(1*Lc+.05*La)*X';
        B2=X*(1*Lm+.05*Ls)*X';
        [V,D]=eig(A2,B2);

        [sorted,index] = sort(diag(D),'descend');
        W=V;
        for i=1:size(index)
            W(:,i)=V(:,index(i));
        end
        V=W;
        V=V(:,1:dim);

        z=V'*X;

        X2=sum(z.^2,1);
        correct=0;
        for i=1:size(yt)
            test=V'*XT(:,i);
            x2=sum(test.^2);
            distance=repmat(x2,1,N)+X2-2*test'*z;
            [sorted,index] = sort(distance,'ascend');
            if yt(i)==y(index(1))
                correct=correct+1;
            end
        end
        precision=correct/size(yt,1);
        myprecision=myprecision+precision;
        disp([' pricision:' num2str(precision)  '  round:' num2str(round)]);



        %     % Zhang's method
        %     Q=ones(N)/N^2;
        %     for i=1:N
        %         for j=1:N
        %             if M(i,j)==1
        %                 Q(i,j)=Q(i,j)-20/nm;
        %             end
        %             if C(i,j)==1
        %                 Q(i,j)=Q(i,j)+1/nc;
        %             end
        %         end
        %     end
        %     Dq = sum(Q(:,:),2);
        %     Lq = diag(Dq)-Q;
        %
        %     A3=X*Lq*X';
        %     [V,D]=eig(A3);
        %     [sorted,index] = sort(diag(D),'descend');
        %     W=V;
        %     for i=1:size(index)
        %         W(:,i)=V(:,index(i));
        %     end
        %     V=W;
        %     V=V(:,1:dim);
        %
        %     z=V'*X;
        %
        %     X2=sum(z.^2,1);
        %     correct=0;
        %     for i=1:size(yt)
        %         test=V'*XT(:,i);
        %         x2=sum(test.^2);
        %         distance=repmat(x2,1,N)+X2-2*test'*z;
        %         [sorted,index] = sort(distance,'ascend');
        %         if yt(i)==y(index(1))
        %             correct=correct+1;
        %         end
        %     end
        %     precision=correct/size(yt,1);
        %     zhangprecision=zhangprecision+precision;
    end




    % % Baseline method
    % X=x';%训练集
    % XT=xt';%测试集
    %
    % X2=sum(X.^2,1);
    % correct=0;
    % for i=1:size(yt)
    %     test=XT(:,i);
    %     x2=sum(test.^2);
    %     distance=repmat(x2,1,N)+X2-2*test'*X;
    %     [sorted,index] = sort(distance,'ascend');
    %     if yt(i)==y(index(1))
    %         correct=correct+1;
    %     end
    % end
    % baseprecision=correct/size(yt,1);
    %
    %
    % % pca method
    % X=x';%训练集
    % XT=xt';%测试集
    %
    % options1.ReducedDim = dim;
    % [eigvector, eigvalue] = PCA(X', options1);
    % Xtrain = X'*eigvector;
    % Xtest=XT'*eigvector;
    % Xtrain=Xtrain';
    % Xtest=Xtest';
    %
    % X2=sum(Xtrain.^2,1);
    % correct=0;
    % for i=1:size(yt)
    %     test=Xtest(:,i);
    %     x2=sum(test.^2);
    %     distance=repmat(x2,1,N)+X2-2*test'*Xtrain;
    %     [sorted,index] = sort(distance,'ascend');
    %     if yt(i)==y(index(1))
    %         correct=correct+1;
    %     end
    % end
    % pcaprecision=correct/size(yt,1);
    %
    %
    % % lda method
    % X=x';%训练集
    % XT=xt';%测试集
    %
    % options2.PCARatio = 0.98;
    % [eigvector, eigvalue] = PCA(X', options2);
    % X = X'*eigvector;
    % XT=XT'*eigvector;
    % X=X';
    % XT=XT';
    %
    % [eigvector, eigvalue, Y] = LDA(X', y);
    % Xtest=XT'*eigvector;
    % Xtrain=Y';
    % Xtest=Xtest';
    %
    % X2=sum(Xtrain.^2,1);
    % correct=0;
    % for i=1:size(yt)
    %     test=Xtest(:,i);
    %     x2=sum(test.^2);
    %     distance=repmat(x2,1,N)+X2-2*test'*Xtrain;
    %     [sorted,index] = sort(distance,'ascend');
    %     if yt(i)==y(index(1))
    %         correct=correct+1;
    %     end
    % end
    % ldaprecision=correct/size(yt,1);

    disp(['Number of K: ' num2str(K)]);
    disp(['Number of target dimensionality: ' num2str(dim)]);
    disp(['Number of constraints: ' num2str(NOC)]);
    disp(['Precision of lgssdr: ' num2str(myprecision/run)]);
    % disp(['Precision of zhang method: ' num2str(zhangprecision/run)]);
    % disp(['Precision of baseline method: ' num2str(baseprecision)]);
    % disp(['Precision of pca method: ' num2str(pcaprecision)]);
    % disp(['Precision of lda method: ' num2str(ldaprecision)]);
end