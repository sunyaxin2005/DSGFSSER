function [ p  ] = locate_face_reco( XL, L, XUL, XTest, LTest)
%LOCATE_FACE_RECO Summary of this function goes here
%   Detailed explanation goes here
[XL, XUL, XTest] = semi_dimension_reduction_main(XL, L, XUL, XTest, 'LDA', 100);
UL = unique(L);
NTest = size(XTest, 1);
pknn = [];
classess = [];
classess1 = LTest; 
for ii = 1 : NTest
    x = XTest(ii, :);
    [ XLL, LL] = find_nearst_classe( XL, L, x, 30);
    [XLL, ~, x] = semi_dimension_reduction_main(XLL, LL, XUL, x, 'LDA', 100);
    [ classes ] = my1nn( XLL, LL, x);
    classess = [classess;classes];
%     sign = 0;
%     ULL = unique(LL);
%     for jj = 1 : 10
%         if(LTest(ii) == ULL(jj))
%            classess1(ii) =  LTest(ii);
%            sign = 1;
%         end
%     end
%     if sign ~= 1
%         classess1(ii) =  0;
%     end
    %p = length(LTest(classess == LTest(1:ii)))/length(classess);
end
p = length(LTest(classess == LTest))/length(classess)
pknn = [pknn;p];
end

function [ XLL, LL, IX1] = find_nearst_classe( XL, L, x, k)
UL = unique(L);
diss = [];
for ii = 1 : length(UL)
    Xthis = XL(L == UL(ii), :);
    dis = pdist2(x, Xthis);
    diss = [diss;min(dis)];
end
[B,IX] = sort(diss);
XLL = [];
LL = [];
for ii = 1 : k
    XLL = [XLL;XL(L == UL(IX(ii)), :)];
    LL = [LL;L(L == UL(IX(ii)))];
end
end