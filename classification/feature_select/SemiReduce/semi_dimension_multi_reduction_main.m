function [ eigvector, Ws ] = semi_dimension_multi_reduction_main( XL, L, XUL, XTest, ALG_TYPE, para )
%SEMI_DIMENSION_MULTI_REDUCTION 此处显示有关此函数的摘要
%   此处显示详细说明
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
    para.k = K1;
    para.K = K1;
end
num = 1;

switch(ALG_TYPE)
     case 'LRGPSSDRW'
        options.PCARatio = 1.0;
        X = [XL;XUL];
        [eigvector, eigvalue] = PCA(X, options);
        XL = XL * eigvector;
        if ~isempty(XUL)
            XUL = XUL * eigvector;
        end
        %eigvector = 1;
        [ Ws ] = LRGPSSDRW_Multi_subspace( XL, L, XUL,  num);
    case 'LRGNEWLLE'
        options.PCARatio = 1;
        X = [XL;XUL];
        [eigvector, eigvalue] = PCA(X, options);
        XL = XL * eigvector;
        if ~isempty(XUL)
            XUL = XUL * eigvector;
        end
        %eigvector = 1;
        [ Ws ] = Multi_LRGNEWLLE( XL, L, XUL, num, para  );
        %[Ws, ~] = LRGNEWLLE_Multi_FS(XL, L, XUL, para);
    case {'LRGPSSDR_GCMCTRACE', 'LRGPSSDR_GCMCTRACE1', 'LRGPSSDR_GCMCTRACE2'}
        options.PCARatio = 1.0;
        X = [XL;XUL];
        [eigvector, eigvalue] = PCA(X, options);
        XL = XL * eigvector;
        if ~isempty(XUL)
            XUL = XUL * eigvector;
        end
        [ Ws ] = Multi_LRGPSSDR_GCMCTRACE( XL, L, XUL, num );  
    case 'EPPFGEDR_GCMC'
         options.PCARatio = 1.0;
        X = [XL;XUL];
        [eigvector, eigvalue] = PCA(X, options);
        XL = XL * eigvector;
        if ~isempty(XUL)
            XUL = XUL * eigvector;
        end
       [ Ws  ] = Multi_EPPFGEDR_GCMC(  XL, L, XUL, num);      
    case 'LRGPSSDR'
          options.PCARatio = 1.0;
        X = [XL;XUL];
        [eigvector, eigvalue] = PCA(X, options);
        XL = XL * eigvector;
        if ~isempty(XUL)
            XUL = XUL * eigvector;
        end
       [ Ws  ] = Multi_LRGPSSDR(  XL, L, XUL, num);          
end

end

