function [eigvector, W] = semi_dimension_reduction_main( XL, L, XUL, XTest, ALG_TYPE, para)
%SEMI_DIMENSION_REDUCTION_MAIN Summary of this function goes here
%   Detailed explanation goes here
UL = unique(L);
K1 = 10000;
for z = 1 : length(UL)
    if ( K1 > length(L(L == UL(z))))
         K1 = length(L(L == UL(z)));
    end
end
 
if nargin < 6
    para.a = 0.1;
    para.b = 0.1;
    para.k = 5;
end

        
switch(ALG_TYPE)
    case 'SSLFDA'
        options.PCARatio = 1.0;
        X = [XL;XUL];
        [eigvector, eigvalue] = PCA(X, options);
        XL = XL * eigvector;
        XUL = XUL * eigvector;
        [ W ] = SSLFDA( XL, L, XUL );
    case 'SDA'
        eigvector = 1;
        [W] = SDA_MAIN(XL, L, XUL, para);
    case 'LRGPSSDA'%目标函数不能奇异
        options.PCARatio = 1.0;
        X = [XL;XUL];
        [eigvector, eigvalue] = PCA(X, options);
        XL = XL * eigvector;
        if ~isempty(XUL)
            XUL = XUL * eigvector;
        end
        [ W ] = LRGPSSDA( XL, L, XUL );
    case 'LRGPSSDR'%目标函数不能奇异
        options.PCARatio = 1.0;
        X = [XL;XUL];
        [eigvector, eigvalue] = PCA(X, options);
        XL = XL * eigvector;
        if ~isempty(XUL)
            XUL = XUL * eigvector;
        end
        [ W ] = LRGPSSDR( XL, L, XUL );
    case 'LRGPSSDR_ITER'%%无用
         options.PCARatio = 1.0;
         X = [XL;XUL];
         [eigvector, eigvalue] = PCA(X, options);
         XL = XL * eigvector;
         if ~isempty(XUL)
            XUL = XUL * eigvector;
         end
         [ W ] = LRGPSSDR_ITER( XL, L, XUL );
    case 'LFDA'%目标函数不能奇异
        options.PCARatio = 1.0;
        X = [XL;XUL];
        [eigvector, eigvalue] = PCA(X, options);
        XL = XL * eigvector;
        
        [W,Z]=LFDA(XL',L, size(XL, 2),  'plain' , para.k);
    case 'PCA'
        options.PCARatio = 1.0;
        X = [XL;XUL];
        eigvector = 1;
        [W, D] = PCA(X, options);
%         [sorted,index] = sort(diag(D),'descend');
%         V = W;
%         for i=1:size(index)
%             W(:,i) = V(:,index(i))./(sum(V(:,index(i))));
%         end
    case 'LDA'%目标函数不能奇异，此函数在内部进行了PCA变换
        options = [];
        options.Fisherface = 1;
        eigvector = 1;
        [W, eigvalue] = LDA(L, options, XL);
    case 'LSDA'
        options = [];
        options.k = para.k;
        eigvector = 1;
        [W, eigvalue] = LSDA(L, options, XL);
    case 'NPE'
        options = [];
        options.k = para.k;
        options.NeighborMode = 'Supervised';
        options.gnd = L;
        eigvector = 1;
        options.ReducedDim = 1000;
        [W, eigvalue] = NPE(options, XL);
    case 'LPP'
        eigvector = 1;
        %X = [XL;XUL];
%       options = [];
     % options.Metric = 'Euclidean';
      options.NeighborMode = 'Supervised';
      options.gnd = L;
      options.bLDA = 1;
      options.ReducedDim = 1000;
      options.k = para.k;
      W = constructW(XL,options);      
      options.PCARatio = 1;
      [W, eigvalue] = LPP(W, options, XL);
    case 'MMP'
        options = [];
        options.WOptions = [];
        eigvector = 1;
        [W] = MMP_MAIN(XL, L, XUL);
    case 'locateLRGPSSDR'
        options.PCARatio = 1.0;
        X = [XL;XUL];
        [eigvector, eigvalue] = PCA(X, options);
        XL1 = XL * eigvector;
        XUL1 = XUL * eigvector;
        [ W ] = LRLFSDR( XL1, L, XUL1);
    case 'LRGAAFF'
        X = [XL;XUL];
         options.PCARatio = 1.0;
        [eigvector, eigvalue] = PCA(X, options);
        XL = XL * eigvector;
        XUL = XUL * eigvector;
        [ W ] = LRGAAFF( XL, L, XUL, 3);
    case 'LRGAAFFSRC'
        X = [XL;XUL];
         options.PCARatio = 1.0;
        [eigvector, eigvalue] = PCA(X, options);
        XL1 = XL * eigvector;
        XUL1 = XUL * eigvector;
        [ W ] = LRGAAFFSRC( XL1, L, XUL1, 3);
    case 'LMNN_SWST'
        options.PCARatio = 1.0;
        X = [XL;XUL];
        [eigvector, eigvalue] = PCA(X, options);
        XL = XL * eigvector;
        [ W ] = LMNN_SWST( XL, L );
    case 'LRGPSPSSDA'
        options.PCARatio = 1.0;
        X = [XL;XUL];
        [eigvector, eigvalue] = PCA(X, options);
        XL = XL * eigvector;
        XUL = XUL * eigvector;
        [ W ] = LRGPSPSSDA( XL, L, XUL );
    case 'LRGPSSDR_SPARSE'
        options.PCARatio = 1.0;
        X = [XL;XUL];
        [eigvector, eigvalue] = PCA(X, options);
        XL = XL * eigvector;
        XUL = XUL * eigvector;
        [ W ] = LRGPSSDR_SPARSE(  XL, L, XUL );
    case 'LRGPSSDRCOS'
        options.PCARatio = 1.0;
        X = [XL;XUL];
        [eigvector, eigvalue] = PCA(X, options);
        XL = XL * eigvector;
        XUL = XUL * eigvector;
        [ W ] = LRGPSSDRCOS( XL, L, XUL );
    case 'LRGPSSDRW'
        options.PCARatio = 1.0;
        X = [XL;XUL];
        [eigvector, eigvalue] = PCA(X, options);
        XL = XL * eigvector;
        if ~isempty(XUL)
            XUL = XUL * eigvector;
        end
        [ W ] = LRGPSSDRW( XL, L, XUL );
    case 'LRGNEWLLE'
        options.PCARatio = 1.0;
        X = [XL;XUL];
        [eigvector, eigvalue] = PCA(X, options);
        XL = XL * eigvector;
        if ~isempty(XUL)
            XUL = XUL * eigvector;
        end
        [ W ] = LRGNEWLLE( XL, L, XUL, para  );
    case 'L21FLDA'
        [ W ] = L21FLDA( XL, L );
        eigvector = 1;
    case 'DSNPE'
        options.PCARatio = 1.0;
        X = [XL;XUL];
        [eigvector, eigvalue] = PCA(X, options);
        eigvector =  eigvector(:, 1: (size(XL, 1)- length(unique(L))));
        XL = XL * eigvector;
        XUL = XUL * eigvector;
%eigvector = 1;
        [ W ] = DSNPE(  XL, L, XUL);
    case 'SPP_train'
        [eigvector, W] = SPP_train( XL, L );
    case 'DSPP_train'
        [eigvector, W] = DSPP_train( XL, L );
    case 'OPSRDR'
        [ eigvector, W ] = OPSRDR( XL, L );
    case 'DimR_SR_New'
         X = [XL;XUL;XTest];
         options.PCARatio = 1.0;
         [eigvector, eigvalue] = PCA(X, options);
        [W]=DimR_SR_New(X',size(eigvector, 2),0.005,2,15);
    case 'LRGPSSDR_GCMC'
        options.PCARatio = 1.0;
        X = [XL;XUL];
        [eigvector, eigvalue] = PCA(X, options);
        XL = XL * eigvector;
        if ~isempty(XUL)
            XUL = XUL * eigvector;
        end
        [ W ] = LRGPSSDR_GCMC( XL, L, XUL );
    case 'LRGPSSDR_NEW_LLE_AUTHOR'
        options.PCARatio = 1.0;
        X = [XL;XUL];
        [eigvector, eigvalue] = PCA(X, options);
        XL = XL * eigvector;
        if ~isempty(XUL)
            XUL = XUL * eigvector;
        end
       [ W  ] = LRGPSSDR_NEW_LLE_AUTHOR(  XL, L, XUL, para );
     case 'LRGPSSDR_NEW_LLE_AUTHOR1'
        options.PCARatio = 1.0;
        X = [XL;XUL];
        [eigvector, eigvalue] = PCA(X, options);
        XL = XL * eigvector;
        if ~isempty(XUL)
            XUL = XUL * eigvector;
        end
       [ W  ] = LRGPSSDR_NEW_LLE_AUTHOR1(  XL, L, XUL, para);      
end

end

function [W] = SDA_MAIN(XL, L, XUL, para)
fea = [XL;XUL];
semiSplit = zeros(1, size(fea, 1));
semiSplit(1 : length(L)) = 1;
semiSplit = (semiSplit > 0);
gnd = zeros(1,  size(fea, 1));
gnd(1 : length(L)) = L;
 
options = [];  
options.k = para.k;

[V, D] = SDA(gnd,fea,semiSplit,options);
[sorted,index] = sort(diag(D),'descend');
W=V;
for i=1:size(index)
    W(:,i)=V(:,index(i));
end
end

function [W] = MMP_MAIN(XL, L, XUL)
fea = [XL;XUL];

semiSplit = zeros(1, size(fea, 1));
semiSplit(1 : length(L)) = 1;
semiSplit = (semiSplit > 0);
gnd = zeros(1,  size(fea, 1));
gnd(1 : length(L)) = L;
[W] = locate_affinity_SL(XL', L, XUL', 3);
 
options = [];  
options.W = W;  
options.ReguBeta = 1; 
options.ReguAlpha = 0.1; 
[V, D] = MMP(L,XL,XUL,options);
[sorted,index] = sort(diag(D),'descend');
W=V;
for i=1:size(index)
    W(:,i)=V(:,index(i));
end
end

function [W] = locate_affinity_SL(XL, L, XUL, K)
%第一步：求局部映射矩阵A
X = [XL, XUL];
X = X';
distance = pdist2(X, X);
[sorted,index] = sort(distance);
neighborhood = index(2:(1+K),:);
sita = sorted(K+1, :);
sita = sita' *  sita;
%distance1 = zeros(size(distance));
% for ii = 1 : size(distance, 2)
%     distance1(neighborhood(:, ii), ii) = distance(neighborhood(:, ii), ii);
%     distance1(ii, neighborhood(:, ii)) = distance(neighborhood(:, ii), ii);
% end
A = exp(-(distance.^2)./sita);
diagA = diag(ones(size(distance, 2), 1));
A = A - diagA;
%A(distance1 == 0) = 0;
X = X';
[D,N] = size(X);
[LD, LN] = size(XL);
%第二步：求W,
W = zeros(N, N);
UL = unique(L);
nClass = length(UL);

%赋值有label的
for ii = 1 : nClass
    W(L == UL(ii), L == UL(ii)) =  A(L == UL(ii), L == UL(ii))/length(L(L == UL(ii)));
end% %赋值没有label的
for ii = 1 : N
    W(neighborhood(:, ii), ii) = W(neighborhood(:, ii), ii) + A(neighborhood(:, ii), ii)/N;
    W(ii, neighborhood(:, ii)) = W(ii, neighborhood(:, ii)) + A(ii, neighborhood(:, ii))/N;
end
% for ii = 1 : N
%     W(neighborhood(:, ii), ii) = A(neighborhood(:, ii), ii)/N;
%     W(ii, neighborhood(:, ii)) = A(ii, neighborhood(:, ii))/N;
% end
%其它的为零
end