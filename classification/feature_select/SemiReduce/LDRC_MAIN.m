function [ eigvector, W, H ] = LDRC_MAIN( XL, L, XUL, XTest, ALG_TYPE, para )
%LDRC_MAIN 此处显示有关此函数的摘要
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
end

        
switch(ALG_TYPE)
    case 'WLDRC1'
        options.PCARatio = 1.0;
        X = [XL;XUL];
        [eigvector, eigvalue] = PCA(X, options);
        XL = XL * eigvector;
        if ~isempty(XUL)
            XUL = XUL * eigvector;
        end
        [ W, H ] = WLDRC1( XL, L, XUL );
    case 'SWLDRC1'
        options.PCARatio = 1.0;
        X = [XL;XUL];
        [eigvector, eigvalue] = PCA(X, options);
        XL = XL * eigvector;
        if ~isempty(XUL)
            XUL = XUL * eigvector;
        end
        [ W, H ] = SWLDRC1( XL, L, XUL );  
    case 'LDRC1'
        options.PCARatio = 1.0;
        X = [XL;XUL];
        [eigvector, eigvalue] = PCA(X, options);
        XL = XL * eigvector;
        if ~isempty(XUL)
            XUL = XUL * eigvector;
        end
        [ W, H ] = LDRC1( XL, L, XUL );   
    case 'LDRCLRGPSSDE1'
        options.PCARatio = 1.0;
        X = [XL;XUL];
        [eigvector, eigvalue] = PCA(X, options);
        XL = XL * eigvector;
        if ~isempty(XUL)
            XUL = XUL * eigvector;
        end
        [ W, H ] = LDRCLRGPSSDE1( XL, L, XUL ); 
end

end

