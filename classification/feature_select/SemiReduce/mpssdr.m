function [ output_args ] = mpssdr( input_args )
%MPSSDR Summary of this function goes here
%   Detailed explanation goes here
for dimT=100:10:100 %目标维数取值范围
    for tr=10:10:50 %训练样本数取值范围
        for knn=3:1:3 %K取值范围
            for NOC=100:100:100 %NOC取值范围
                
                %初始化参数
                mpssdrprecision1=0;
                mpssdrCMGprecision1=0;
                mpssdrCMLprecision1=0;
                mpssdrCMprecision1=0;
                mpssdrLprecision1=0;
                mpssdrGprecision1=0;
                npssdrprecision1=0;
                lgssdrprecision1=0;
                ssdrprecision1=0;
                clppprecision1=0;
                baseprecision1=0;
                pcaprecision1=0;
                run=5;%实验次数
                nmNumber=0;%如果为0，则随机生成约束
                
                Kname=['result/YaleB_32x32_date=121225_tr=',num2str(tr)];
                Kname=[Kname,'_k='];
                Kname=[Kname,num2str(knn)];
                Kname=[Kname,'_dim='];
                Kname=[Kname,num2str(dimT)];
                Kname=[Kname,'_NOC='];
                Kname=[Kname,num2str(NOC)];
                
                Kname=[Kname,'.txt'];
                fid = fopen(Kname, 'w');
                                
                for indexARG=1:20
                    FileName=['data\YaleB_32x32_',num2str(tr)];
                    FileName=[FileName,'Train_'];
                    FileName=[FileName,num2str(indexARG)];
                    FileName=[FileName,'.mat'];
                    load(FileName);
                    
                    X=double(x);%训练集
                    XT=double(xt);%测试集
                    N=length(y);%训练数据个数
                    clear x xt;
                    
                    dim=dimT;%目标维数
                    K=knn;%LLE最近邻接点个数(不包括自身)
                                        
                    X=X';
                    XT=XT';
                    baseprecision=0;
                    pcaprecision=0;
                    % Baseline
                    X2=sum(X.^2,1);
                    correct=0;
                    for i=1:length(yt)
                        test=XT(:,i);
                        x2=sum(test.^2);
                        distance=repmat(x2,1,N)+X2-2*test'*X;
                        [sorted,index] = sort(distance,'ascend');
                        if yt(i)==y(index(1))
                            correct=correct+1;
                        end
                    end
                    baseprecision=baseprecision+correct/length(yt);
                    baseprecision1=baseprecision1+baseprecision;
                    clear distance sorted index test x2 X2  i correct;
                    
                    % PCA
                    options1.ReducedDim = dim;
                    [eigvector, eigvalue] = PCA(X', options1);
                    
                    Xtrain = X'*eigvector;
                    Xtest=XT'*eigvector;
                    Xtrain=Xtrain';
                    Xtest=Xtest';
                    
                    X2=sum(Xtrain.^2,1);
                    correct=0;
                    for i=1:length(yt)
                        test=Xtest(:,i);
                        x2=sum(test.^2);
                        distance=repmat(x2,1,size(Xtrain,2))+X2-2*test'*Xtrain;
                        [sorted,index] = sort(distance,'ascend');
                        if yt(i)==y(index(1))
                            correct=correct+1;
                        end
                    end
                    pcaprecision=pcaprecision+correct/length(yt);
                    pcaprecision1=pcaprecision1+pcaprecision;
                    clear Xtest Xtrain distance index sorted x2 test eigvector  eigvalue;
                    X=X';
                    XT=XT';
                                        
                    options.PCARatio = 0.98;
                    [eigvector, eigvalue] = PCA(X, options);
                    X = X*eigvector;
                    XT=XT*eigvector;
                    X=X';
                    XT=XT';
                    
                    if K~=0
                        [Dt,Nt] = size(X);
                        % STEP1: COMPUTE PAIRWISE DISTANCES & FIND NEIGHBORS
                        X2 = sum(X.^2,1);
                        distance = repmat(X2,Nt,1)+repmat(X2',1,Nt)-2*X'*X;
                        clear X2;
                        [sorted,index] = sort(distance);
                        neighborhood = index(2:(1+K),:);
                        
                        % STEP2: SOLVE FOR RECONSTRUCTION WEIGHTS
                        if(K>Dt)
                            tol=1e-3; % regularlizer in case constrained fits are ill conditioned
                        else
                            tol=1e-12;
                        end
                        Wt = zeros(K,Nt);
                        for ii=1:Nt
                            z = X(:,neighborhood(:,ii))-repmat(X(:,ii),1,K); % shift ith pt to origin
                            Ct = z'*z;                                        % local covariance
                            Ct = Ct + eye(K,K)*tol*trace(Ct);                   % regularlization (K>Dt)
                            Wt(:,ii) = Ct\ones(K,1);                           % solve Cw=1
                            Wt(:,ii) = Wt(:,ii)/sum(Wt(:,ii));                  % enforce sum(w)=1
                        end;
                        
                        SW=sparse(zeros(N));
                        for i=1:N
                            for j=2:1+K
                                SW(i,index(j,i))=Wt(j-1,i);
                            end
                        end
                        clear distance sorted index;
                        clear SW;%没有用的时候暂时不用
                        % STEP 3: COMPUTE EMBEDDING FROM EIGENVECTS OF COST MATRIX Mt=(I-Wt)'(I-Wt)
                        Mt = sparse(1:Nt,1:Nt,ones(1,Nt),Nt,Nt,4*K*Nt);
                        for ii=1:Nt
                            w = Wt(:,ii);
                            jj = neighborhood(:,ii);
                            Mt(ii,jj) = Mt(ii,jj) - w';
                            Mt(jj,ii) = Mt(jj,ii) - w;
                            Mt(jj,jj) = Mt(jj,jj) + w*w';
                        end;
                        clear w neighborhood;
                        % options.k = K;
                        % [fTauDg] = IsoP(options, X');
                    else
                        Mt=0;
                    end                
                    
                    mpssdrprecision=0;
                    mpssdrCMGprecision=0;
                    mpssdrCMLprecision=0;
                    mpssdrCMprecision=0;
                    mpssdrLprecision=0;
                    mpssdrGprecision=0;
                    npssdrprecision=0;
                    lgssdrprecision=0;
                    ssdrprecision=0;
                    clppprecision=0;
                    
                    for round=1:run
                        M=(zeros(N,N));%must-link矩阵
                        C=(zeros(N,N));%cannot-link矩阵
                        nm=0;%must-link点对个数
                        nc=0;%cannot-link点对个数
                        S=(eye(N,N));%邻接矩阵
                        A=ones(N,N)-eye(N,N);%非邻接矩阵
                                                
                        if nmNumber==0  %随机生成约束
                            for i=1:NOC
                                part1=randint(1,1,[1,N]);
                                part2=randint(1,1,[1,N]);
                                while part1==part2 || M(part1,part2)==1 || C(part1,part2)==1
                                    part2=randint(1,1,[1,N]);
                                end
                                if  y(part1)==y(part2)
                                    M(part1,part2)=1;
                                    M(part2,part1)=1;
                                    nm=nm+1;
                                end
                                if  y(part1)~=y(part2)
                                    C(part1,part2)=1;
                                    C(part2,part1)=1;
                                    nc=nc+1;
                                end
                            end
                        else
                            while(nm<nmNumber||nc<NOC-nmNumber)
                                part1=randint(1,1,[1,N]);
                                part2=randint(1,1,[1,N]);
                                while part1==part2 || M(part1,part2)==1 || C(part1,part2)==1
                                    part2=randint(1,1,[1,N]);
                                end
                                if nm<nmNumber && y(part1)==y(part2)
                                    M(part1,part2)=1;
                                    M(part2,part1)=1;
                                    nm=nm+1;
                                end
                                if nc<NOC-nmNumber && y(part1)~=y(part2)
                                    C(part1,part2)=1;
                                    C(part2,part1)=1;
                                    nc=nc+1;
                                end
                            end
                        end
                        
                        X2 = sum(X.^2,1);
                        distance1 = repmat(X2,N,1)+repmat(X2',1,N)-2*X'*X;
                        [sorted1,index1] = sort(distance1);
                        
                        K1=K;
                        for i=1:N
                            for j=2:1+K1
                                if M(i,index1(j,i))==0 && C(i,index1(j,i))==0
                                    S(i,index1(j,i))=1;
                                    S(index1(j,i),i)=1;
                                    A(i,index1(j,i))=0;
                                    A(index1(j,i),i)=0;
                                end
                            end
                        end
                        clear distance1 sorted1 index1;
                        
                        %                         Ds = sum(S(:,:),2);
                        %                         Ls = diag(Ds)-S;
                        clear S;
                        Da = sum(A(:,:),2);
                        La = diag(Da)-A;
                        clear A;
                        Dm = sum(M(:,:),2);
                        Lm = diag(Dm)-M;
                        Dc = sum(C(:,:),2);
                        Lc = diag(Dc)-C;
                        
                        %                         H1=1/(2*nc);
                        %                         H2=1/(1*N*(N-K));
                        %                         H3=100/(2*nm);
                        %                         H4=1;
                        %                         H5=H1+H2+H3+H4;
                        %                         H6=H1+H2+H4;
                        %                         H1=H1/H5;H11=H1/H6;
                        %                         H2=H2/H5;H22=H2/H6;
                        %                         H3=H3/H5;
                        %                         H4=H4/H5;H44=H4/H6;
                        
                        % MPSSDR
                        if nm~=0
                            %                             A=X*(H1*Lc+H2*La-H3*Lm-H4*Mt)*X';
                            A=X*(1*Lc/(1*nc)+1*La/(1*N*(N-K))-20*Lm/(1*nm)-1000*Mt/(1*N))*X';
                        else
                            %                             A=X*(H11*Lc+H22*La-H44*Mt)*X';
                            A=X*(1*Lc/(1*nc)+1*La/(1*N*(N-K))-1000*Mt/(1*N))*X';                           
                        end
                        [V,D]=eig(A);
                        
                        [sorted,index] = sort(diag(D),'descend');
                        W=V;
                        for i=1:length(index)
                            W(:,i)=V(:,index(i));
                        end
                        V=W;
                        V=V(:,1:dim);
                        
                        z=V'*X;
                        
                        X2=sum(z.^2,1);
                        correct=0;
                        for i=1:length(yt)
                            test=V'*XT(:,i);
                            x2=sum(test.^2);
                            distance=repmat(x2,1,N)+X2-2*test'*z;
                            [sorted,index] = sort(distance,'ascend');
                            if yt(i)==y(index(1))
                                correct=correct+1;
                            end
                        end
                        precision=correct/length(yt);
                        mpssdrprecision=mpssdrprecision+precision;
                        clear test x2 diatance sorted index;
                        
                        % MPSSDR-CMG 约束信息+全局结构信息
                        if nm~=0
                            %                             A=X*(H1*Lc+H2*La-H3*Lm-H4*Mt)*X';
                            A=X*(1*Lc/(1*nc)+1*La/(1*N*(N-K))-20*Lm/(1*nm)-0*Mt/(1*N))*X';
                        else
                            %                             A=X*(H11*Lc+H22*La-H44*Mt)*X';
                            A=X*(1*Lc/(1*nc)+1*La/(1*N*(N-K))-0*Mt/(1*N))*X';                       
                        end
                        [V,D]=eig(A);
                        
                        [sorted,index] = sort(diag(D),'descend');
                        W=V;
                        for i=1:length(index)
                            W(:,i)=V(:,index(i));
                        end
                        V=W;
                        V=V(:,1:dim);
                        
                        z=V'*X;
                        
                        X2=sum(z.^2,1);
                        correct=0;
                        for i=1:length(yt)
                            test=V'*XT(:,i);
                            x2=sum(test.^2);
                            distance=repmat(x2,1,N)+X2-2*test'*z;
                            [sorted,index] = sort(distance,'ascend');
                            if yt(i)==y(index(1))
                                correct=correct+1;
                            end
                        end
                        precision=correct/length(yt);
                        mpssdrCMGprecision=mpssdrCMGprecision+precision;
                        clear test x2 diatance sorted index;
                        
                        % MPSSDR-CML 约束信息+局部结构信息
                        if nm~=0
                            %                             A=X*(H1*Lc+H2*La-H3*Lm-H4*Mt)*X';
                            A=X*(1*Lc/(1*nc)+0*La/(1*N*(N-K))-20*Lm/(1*nm)-1000*Mt/(1*N))*X';
                        else
                            %                             A=X*(H11*Lc+H22*La-H44*Mt)*X';
                            A=X*(1*Lc/(1*nc)+0*La/(1*N*(N-K))-1000*Mt/(1*N))*X';                          
                        end
                        [V,D]=eig(A);
                        
                        [sorted,index] = sort(diag(D),'descend');
                        W=V;
                        for i=1:length(index)
                            W(:,i)=V(:,index(i));
                        end
                        V=W;
                        V=V(:,1:dim);
                        
                        z=V'*X;
                        
                        X2=sum(z.^2,1);
                        correct=0;
                        for i=1:length(yt)
                            test=V'*XT(:,i);
                            x2=sum(test.^2);
                            distance=repmat(x2,1,N)+X2-2*test'*z;
                            [sorted,index] = sort(distance,'ascend');
                            if yt(i)==y(index(1))
                                correct=correct+1;
                            end
                        end
                        precision=correct/length(yt);
                        mpssdrCMLprecision=mpssdrCMLprecision+precision;
                        clear test x2 diatance sorted index;
                        
                        % MPSSDR-CM 约束信息
                        if nm~=0
                            %                             A=X*(H1*Lc+H2*La-H3*Lm-H4*Mt)*X';
                            A=X*(1*Lc/(1*nc)+0*La/(1*N*(N-K))-20*Lm/(1*nm)-0*Mt/(1*N))*X';
                        else
                            %                             A=X*(H11*Lc+H22*La-H44*Mt)*X';
                            A=X*(1*Lc/(1*nc)+0*La/(1*N*(N-K))-0*Mt/(1*N))*X';                          
                        end
                        [V,D]=eig(A);
                        
                        [sorted,index] = sort(diag(D),'descend');
                        W=V;
                        for i=1:length(index)
                            W(:,i)=V(:,index(i));
                        end
                        V=W;
                        V=V(:,1:dim);
                        
                        z=V'*X;
                        
                        X2=sum(z.^2,1);
                        correct=0;
                        for i=1:length(yt)
                            test=V'*XT(:,i);
                            x2=sum(test.^2);
                            distance=repmat(x2,1,N)+X2-2*test'*z;
                            [sorted,index] = sort(distance,'ascend');
                            if yt(i)==y(index(1))
                                correct=correct+1;
                            end
                        end
                        precision=correct/length(yt);
                        mpssdrCMprecision=mpssdrCMprecision+precision;
                        clear test x2 diatance sorted index;
                        
                        % MPSSDR-L 局部结构信息
                        if nm~=0
                            %                             A=X*(H1*Lc+H2*La-H3*Lm-H4*Mt)*X';
                            A=X*(0*Lc/(1*nc)+0*La/(1*N*(N-K))-0*Lm/(1*nm)-1000*Mt/(1*N))*X';
                        else
                            %                             A=X*(H11*Lc+H22*La-H44*Mt)*X';
                            A=X*(0*Lc/(1*nc)+0*La/(1*N*(N-K))-1000*Mt/(1*N))*X';                          
                        end
                        [V,D]=eig(A);
                        
                        [sorted,index] = sort(diag(D),'descend');
                        W=V;
                        for i=1:length(index)
                            W(:,i)=V(:,index(i));
                        end
                        V=W;
                        V=V(:,1:dim);
                        
                        z=V'*X;
                        
                        X2=sum(z.^2,1);
                        correct=0;
                        for i=1:length(yt)
                            test=V'*XT(:,i);
                            x2=sum(test.^2);
                            distance=repmat(x2,1,N)+X2-2*test'*z;
                            [sorted,index] = sort(distance,'ascend');
                            if yt(i)==y(index(1))
                                correct=correct+1;
                            end
                        end
                        precision=correct/length(yt);
                        mpssdrLprecision=mpssdrLprecision+precision;
                        clear test x2 diatance sorted index;
                        
                        % MPSSDR-G 全局结构信息
                        if nm~=0
                            %                             A=X*(H1*Lc+H2*La-H3*Lm-H4*Mt)*X';
                            A=X*(0*Lc/(1*nc)+1*La/(1*N*(N-K))-0*Lm/(1*nm)-0*Mt/(1*N))*X';
                        else
                            %                             A=X*(H11*Lc+H22*La-H44*Mt)*X';
                            A=X*(0*Lc/(1*nc)+1*La/(1*N*(N-K))-0*Mt/(1*N))*X';                          
                        end
                        [V,D]=eig(A);
                        
                        [sorted,index] = sort(diag(D),'descend');
                        W=V;
                        for i=1:length(index)
                            W(:,i)=V(:,index(i));
                        end
                        V=W;
                        V=V(:,1:dim);
                        
                        z=V'*X;
                        
                        X2=sum(z.^2,1);
                        correct=0;
                        for i=1:length(yt)
                            test=V'*XT(:,i);
                            x2=sum(test.^2);
                            distance=repmat(x2,1,N)+X2-2*test'*z;
                            [sorted,index] = sort(distance,'ascend');
                            if yt(i)==y(index(1))
                                correct=correct+1;
                            end
                        end
                        precision=correct/length(yt);
                        mpssdrGprecision=mpssdrGprecision+precision;
                        clear test x2 diatance sorted index;
                        
                        % NPSSDR
                        A=X*(1*Lc)*X';
                        if nm~=0
                            B=X*(1*Lm+.1*Mt)*X';
                        else
                            B=X*(.1*Mt)*X';
                        end
                        [V,D]=eig(A,B);
                        
                        [sorted,index] = sort(diag(D),'descend');
                        W=V;
                        for i=1:length(index)
                            W(:,i)=V(:,index(i));
                        end
                        V=W;
                        V=V(:,1:dim);
                        
                        z=V'*X;
                        
                        X2=sum(z.^2,1);
                        correct=0;
                        for i=1:length(yt)
                            test=V'*XT(:,i);
                            x2=sum(test.^2);
                            distance=repmat(x2,1,N)+X2-2*test'*z;
                            [sorted,index] = sort(distance,'ascend');
                            if yt(i)==y(index(1))
                                correct=correct+1;
                            end
                        end
                        precision=correct/length(yt);
                        npssdrprecision=npssdrprecision+precision;
                        clear  distance  sorted x2 test index;
                                               
                        % SSDR
                        Q=ones(N)/N^2;
                        for i=1:N
                            for j=1:N
                                if M(i,j)==1
                                    Q(i,j)=Q(i,j)-20/nm;
                                end
                                if C(i,j)==1
                                    Q(i,j)=Q(i,j)+1/nc;
                                end
                            end
                        end
                        Dq = sum(Q(:,:),2);
                        Lq = diag(Dq)-Q;
                        clear Q;
                        
                        A=X*Lq*X';
                        [V,D]=eig(A);
                        [sorted,index] = sort(diag(D),'descend');
                        W=V;
                        for i=1:length(index)
                            W(:,i)=V(:,index(i));
                        end
                        V=W;
                        V=V(:,1:dim);
                        
                        z=V'*X;
                        
                        X2=sum(z.^2,1);
                        correct=0;
                        for i=1:length(yt)
                            test=V'*XT(:,i);
                            x2=sum(test.^2);
                            distance=repmat(x2,1,N)+X2-2*test'*z;
                            [sorted,index] = sort(distance,'ascend');
                            if yt(i)==y(index(1))
                                correct=correct+1;
                            end
                        end
                        precision=correct/length(yt);
                        ssdrprecision=ssdrprecision+precision;
                                                
                        % CLPP
                        clear La Lc Lm Ls Lq;
                        S=eye(N);
                        X2 = sum(X.^2,1);
                        distance1 = repmat(X2,N,1)+repmat(X2',1,N)-2*X'*X;
                        [sorted1,index1] = sort(distance1);
                        beta=sum(sum(distance1))/(N*N);
                        for i=1:N
                            for j=2:1+K1
                                S(i,index1(j,i))=exp(-(sorted1(j,i)^2)/(beta^2/1));
                                S(index1(j,i),i)=exp(-(sorted1(j,i)^2)/(beta^2/1));
                            end
                        end
                        
                        WD=(zeros(N));
                        for i=1:N
                            for j=1:N
                                if M(i,j)==1
                                    S(i,j)=1;
                                    S(j,i)=1;
                                    comm=intersect(index1(2:1+K1,i),index1(2:1+K1,j));
                                    if ~isempty(comm)
                                        S(i,comm)=1;
                                        S(comm,i)=1;
                                        S(j,comm)=1;
                                        S(comm,j)=1;
                                    end
                                end
                                if C(i,j)==1
                                    if S(i,j)~=0
                                        S(i,j)=-1;
                                        S(j,i)=-1;
                                        WD(i,j)=-1;
                                        WD(j,i)=-1;
                                    end
                                    comm=intersect(index1(2:1+K1,i),index1(2:1+K1,j));
                                    if ~isempty(comm)
                                        WD(i,j)=-1;
                                        WD(j,i)=-1;
                                        for k=1:length(comm)
                                            if S(i,comm(k))<S(j,comm(k))
                                                S(i,comm(k))=0;
                                                S(comm(k),i)=0;
                                            elseif S(i,comm(k))>S(j,comm(k))
                                                S(j,comm(k))=0;
                                                S(comm(k),j)=0;
                                            end
                                        end
                                    end
                                end
                            end
                        end
                        
                        clear distance1 sorted1 index1;
                        Ds = diag(sum(S(:,:),2));
                        Ls = Ds-S;
                        clear S;
                        Dws = diag(sum(M(:,:),2));
                        Lws = Dws-M;
                        clear M;
                        Dwd = diag(sum(WD(:,:),2));
                        Lwd = Dwd-WD;
                        clear WD;
                        
                        A=X*(Ls+Lws+Lwd)*X';
                        clear Ls Lws Lwd;
                        B=X*(Ds+Dws+Dwd)*X';
                        clear Ds Dws Dwd;
                        [V,D]=eig(A,B);
                        
                        [sorted,index] = sort(diag(D),'ascend');
                        W=V;
                        for i=1:length(index)
                            W(:,i)=V(:,index(i));
                        end
                        V=W;
                        V=V(:,1:dim);
                        
                        z=V'*X;
                        
                        X2=sum(z.^2,1);
                        correct=0;
                        for i=1:length(yt)
                            test=V'*XT(:,i);
                            x2=sum(test.^2);
                            distance=repmat(x2,1,N)+X2-2*test'*z;
                            [sorted,index] = sort(distance,'ascend');
                            if yt(i)==y(index(1))
                                correct=correct+1;
                            end
                        end
                        precision=correct/length(yt);
                        clppprecision=clppprecision+precision;                                                
                    end
                    
                    mpssdrprecision1=mpssdrprecision1+mpssdrprecision/run;
                    mpssdrCMGprecision1=mpssdrCMGprecision1+mpssdrCMGprecision/run;
                    mpssdrCMLprecision1=mpssdrCMLprecision1+mpssdrCMLprecision/run;
                    mpssdrCMprecision1=mpssdrCMprecision1+mpssdrCMprecision/run;
                    mpssdrLprecision1=mpssdrLprecision1+mpssdrLprecision/run;
                    mpssdrGprecision1=mpssdrGprecision1+mpssdrGprecision/run;
                    npssdrprecision1=npssdrprecision1+npssdrprecision/run;
                    lgssdrprecision1=lgssdrprecision1+lgssdrprecision/run;
                    ssdrprecision1=ssdrprecision1+ssdrprecision/run;
                    clppprecision1=clppprecision1+clppprecision/run;
                    disp(['****************************' FileName '****************************']);
                    disp(['Number of K: ' num2str(K)]);
                    disp(['Number of target dimensionality: ' num2str(dim)]);
                    disp(['Number of constraints: ' num2str(NOC)]);
                    disp(['Number of nmNumber: ' num2str(nm)]);
                    disp(['Precision of mpssdr method: ' num2str(mpssdrprecision/run)]);
                    disp(['Precision of mpssdr-CMG method: ' num2str(mpssdrCMGprecision/run)]);
                    disp(['Precision of mpssdr-CML method: ' num2str(mpssdrCMLprecision/run)]);
                    disp(['Precision of mpssdr-CM method: ' num2str(mpssdrCMprecision/run)]);
                    disp(['Precision of mpssdr-L method: ' num2str(mpssdrLprecision/run)]);
                    disp(['Precision of mpssdr-G method: ' num2str(mpssdrGprecision/run)]);
                    disp(['Precision of npssdr method: ' num2str(npssdrprecision/run)]);
                    %     disp(['Precision of lgssdr method: ' num2str(lgssdrprecision/run)]);
                    disp(['Precision of ssdr method: ' num2str(ssdrprecision/run)]);
                    disp(['Precision of clpp method: ' num2str(clppprecision/run)]);
                    disp(['Precision of baseline method: ' num2str(baseprecision)]);
                    disp(['Precision of pca method: ' num2str(pcaprecision)]);
                    disp('****************************');
                    fprintf(fid,'%s \r\n',['************************ ' FileName '***********************']);
                    fprintf(fid,'%s \r\n',['Number of K: ' num2str(K)]);
                    fprintf(fid,'%s \r\n',['Number of target dimensionality: ' num2str(dim)]);
                    fprintf(fid,'%s \r\n',['Number of constraints: ' num2str(NOC)]);
                    fprintf(fid,'%s \r\n',['Number of nmNumber: ' num2str(nm)]);
                    fprintf(fid,'%s \r\n',['Precision of mpssdr method: ' num2str(mpssdrprecision/run)]);
                    fprintf(fid,'%s \r\n',['Precision of mpssdr-CMG method: ' num2str(mpssdrCMGprecision/run)]);
                    fprintf(fid,'%s \r\n',['Precision of mpssdr-CML method: ' num2str(mpssdrCMLprecision/run)]);
                    fprintf(fid,'%s \r\n',['Precision of mpssdr-CM method: ' num2str(mpssdrCMprecision/run)]);
                    fprintf(fid,'%s \r\n',['Precision of mpssdr-L method: ' num2str(mpssdrLprecision/run)]);
                    fprintf(fid,'%s \r\n',['Precision of mpssdr-G method: ' num2str(mpssdrGprecision/run)]);
                    fprintf(fid,'%s \r\n',['Precision of npssdr method: ' num2str(npssdrprecision/run)]);
                    %                     fprintf(fid,'%s \r\n',['Precision of lgssdr method: ' num2str(lgssdrprecision/run)]);
                    fprintf(fid,'%s \r\n',['Precision of ssdr method: ' num2str(ssdrprecision/run)]);
                    fprintf(fid,'%s \r\n',['Precision of clpp method: ' num2str(clppprecision/run)]);
                    fprintf(fid,'%s \r\n',['Precision of baseline method: ' num2str(baseprecision)]);
                    fprintf(fid,'%s \r\n',['Precision of pca method: ' num2str(pcaprecision)]);
                    fprintf(fid,'%s \r\n',['************************************************ ']);                                        
                    clear A B C D Ct Da Dc Dm Dq Dt M Mt V Wt X X2 XT eigvalue eigvector jj index sorted test;
                end
                
                runindex=indexARG;
                disp('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@');
                disp(['Number of K: ' num2str(K)]);
                disp(['Number of target dimensionality: ' num2str(dim)]);
                disp(['Number of constraints: ' num2str(NOC)]);
                disp(['Number of nmNumber: ' num2str(nm)]);
                disp(['Average Precision of mpssdr method: ' num2str(mpssdrprecision1/runindex)]);
                disp(['Average Precision of mpssdr-CMG method: ' num2str(mpssdrCMGprecision1/runindex)]);
                disp(['Average Precision of mpssdr-CML method: ' num2str(mpssdrCMLprecision1/runindex)]);
                disp(['Average Precision of mpssdr-CM method: ' num2str(mpssdrCMprecision1/runindex)]);
                disp(['Average Precision of mpssdr-L method: ' num2str(mpssdrLprecision1/runindex)]);
                disp(['Average Precision of mpssdr-G method: ' num2str(mpssdrGprecision1/runindex)]);
                disp(['Average Precision of npssdr method: ' num2str(npssdrprecision1/runindex)]);
                %     disp(['Average Precision of lgssdr method: ' num2str(lgssdrprecision/run)]);
                disp(['Average Precision of ssdr method: ' num2str(ssdrprecision1/runindex)]);
                disp(['Average Precision of clpp method: ' num2str(clppprecision1/runindex)]);
                disp(['Average Precision of baseline method: ' num2str(baseprecision1/indexARG)]);
                disp(['Average Precision of pca method: ' num2str(pcaprecision1/indexARG)]);
                disp('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@');
                fprintf(fid,'%s \r\n',['@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@']);
                fprintf(fid,'%s \r\n',['Number of K: ' num2str(K)]);
                fprintf(fid,'%s \r\n',['Number of target dimensionality: ' num2str(dim)]);
                fprintf(fid,'%s \r\n',['Number of constraints: ' num2str(NOC)]);
                fprintf(fid,'%s \r\n',['Number of nmNumber: ' num2str(nm)]);
                fprintf(fid,'%s \r\n',['Average Precision of mpssdr method: ' num2str(mpssdrprecision1/runindex)]);
                fprintf(fid,'%s \r\n',['Average Precision of mpssdr-CMG method: ' num2str(mpssdrCMGprecision1/runindex)]);
                fprintf(fid,'%s \r\n',['Average Precision of mpssdr-CML method: ' num2str(mpssdrCMLprecision1/runindex)]);
                fprintf(fid,'%s \r\n',['Average Precision of mpssdr-CM method: ' num2str(mpssdrCMprecision1/runindex)]);
                fprintf(fid,'%s \r\n',['Average Precision of mpssdr-L method: ' num2str(mpssdrLprecision1/runindex)]);
                fprintf(fid,'%s \r\n',['Average Precision of mpssdr-G method: ' num2str(mpssdrGprecision1/runindex)]);
                fprintf(fid,'%s \r\n',['Average Precision of npssdr method: ' num2str(npssdrprecision1/runindex)]);
                %   fprintf(fid,'%s \r\n',['Average Precision of lgssdr method: ' num2str(lgssdrprecision/run)]);
                fprintf(fid,'%s \r\n',['Average Precision of ssdr method: ' num2str(ssdrprecision1/runindex)]);
                fprintf(fid,'%s \r\n',['Average Precision of clpp method: ' num2str(clppprecision1/runindex)]);
                fprintf(fid,'%s \r\n',['Average Precision of baseline method: ' num2str(baseprecision1/indexARG)]);
                fprintf(fid,'%s \r\n',['Average Precision of pca method: ' num2str(pcaprecision1/indexARG)]);
                fprintf(fid,'%s \r\n',['@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@']);
                fclose(fid);
            end
        end
    end
end