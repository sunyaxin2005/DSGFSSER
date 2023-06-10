function [  prs, pstds ] = locate_dimension_reduction( X, L, ALg_type, sel_dims )
%LOCATE_DIMENSION_REDUCTION Summary of this function goes here
%   Detailed explanation goes here
%TEST_SEMI Summary of this function goes here
%   Detailed explanation goes here
%ALg_type = 'LMNN_SWST'
global filename1;
global file_name
ALg_types = {'SSLFDA', 'SDA', 'LRGPSSDR', 'LFDA',  'PCA', 'LDA', 'LSDA', 'NPE', 'LPP', 'MMP', 'locateLRGPSSDR', 'LRGAAFF', 'LRGPSPSSDA', 'LRGPSSDRCOS'};
random_times = 10;
%X = NormalizeFea(X);
%X = NormalizeFea(X);
stdd = std(X);
meandd = mean(X);
stdd = repmat(stdd, size(X, 1), 1);
meandd = repmat(meandd, size(X, 1), 1);
X = (X - meandd)./stdd;
X(isnan(X)) = 0;
X(isinf(X)) = 0;
%X = pdist2(X, X);

p1 = zeros(random_times, length(sel_dims));
p2 = zeros(random_times, length(sel_dims));
p3 = zeros(random_times, length(sel_dims));
p4 = zeros(random_times, length(sel_dims));
p5 = zeros(random_times, length(sel_dims));
p6 = zeros(random_times, length(sel_dims));
p7 = zeros(random_times, length(sel_dims));
p8 = zeros(random_times, length(sel_dims));

for ii = 1 : random_times
    %try
    l =3;
    ul =3;
    us = 7;
    [train, semitrain, test] = get_l_ul_us_inx(L, l, ul);
    %[train, semitrain, test] = get_l_ul_us_inx2(L, 0.5, 0.1);
    XLTrain = X(train, :);
    LTrain = L(train);
    XULTrain = X(semitrain, :);
    XTest = X(test, :);
    LTest = L(test);
 
   
    [eigvector, W] = semi_dimension_reduction_main(XLTrain, LTrain, XULTrain, ALg_type);
    
    for jj = 1 : length(sel_dims)
        str = sprintf('%d_%d', ii, jj);
        filename1 = strcat(file_name, str);
        
        sel_dim = sel_dims(jj);
        if sel_dim > size(W, 2)
            sel_dim = size(W, 2);
        end
        XLTrain1 = XLTrain * eigvector * W(:, 1 : sel_dim);
        XULTrain1 = XULTrain * eigvector * W(:, 1 : sel_dim);
        XTest1 = XTest * eigvector * W(:, 1 : sel_dim);  
        %[ classes ] = LMC1( XTest1', XLTrain1', LTrain, 3);
        [ classes ] = my1nn( XLTrain1, LTrain, XTest1);
        p = length(LTest(classes == LTest))/length(classes);
        p1(ii, jj) = p;
%         [ I1 ] = class2bmp( XLTrain, LTrain );
%         I1 = mat2gray(I1);
%         %figure, imshow(I1);
%         figure, plot(XLTrain(:, 1), XLTrain(:, 2), '*');
%         
%         [ I1 ] = class2bmp( XLTrain1, LTrain );
%         I1 = mat2gray(I1);
%         figure, plot(XLTrain1(:, 1), XLTrain1(:, 2), '*');
        %figure, imshow(I1);
%         model = svmtrain(LTrain,XLTrain1,  '-t polynomial');
%         [classes, accuracy, dec_values] = svmpredict(LTest, XTest1, model);  
%         p = length(LTest(classes == LTest))/length(classes);
%         p2(ii, jj) = p; 
%        
%         [ classes ] = SRC_MAIN( XLTrain1, LTrain, XTest1, 'tradition_src');
%         p = length(LTest(classes == LTest))/length(classes);
%         p2(ii, jj) = p;
% 
%         [ classes ] = SRC_MAIN( XLTrain1, LTrain, XTest1, 'LSDL_SRC');
%         p = length(LTest(classes == LTest))/length(classes);
%         p3(ii, jj) = p;
% %         
%         [ classes ] = SRC_MAIN( XLTrain1, LTrain, XTest1, 'CLSDL_SRC');
%         p = length(LTest(classes == LTest))/length(classes);
%         p4(ii, jj) = p;
% %         
%         [ classes ] = SRC_MAIN( XLTrain1, LTrain, XTest1, 'WSRC');
%         p = length(LTest(classes == LTest))/length(classes);
%         p5(ii, jj) = p;
%         
%         [ classes ] = SRC_MAIN( XLTrain1, LTrain, XTest1, 'LSDL_SRC');
%         p = length(LTest(classes == LTest))/length(classes);
%         p6(ii, jj) = p;
        
%         [ classes ] = SRC_MAIN( XLTrain1, LTrain, XTest1, 'tradition_srcLASSO');
%         p = length(LTest(classes == LTest))/length(classes);
%         p7(ii, jj) = p;
%         
%         [ classes ] = SRC_MAIN( XLTrain1, LTrain, XTest1, 'kcckneighbor_srcLASSO');
%         p = length(LTest(classes == LTest))/length(classes);
%         p8(ii, jj) = p;
        
        
        if sel_dims(jj) > size(W, 2)
            break
        end
    end
%     catch
%         %ii = ii - 1;
%     end
end
p1s = std(p1, 0, 1);
p2s = std(p2, 0, 1);
p3s = std(p3, 0, 1);
p4s = std(p4, 0, 1);
p5s = std(p5, 0, 1);
p6s = std(p6, 0, 1);
p7s = std(p7, 0, 1);
p8s = std(p8, 0, 1);

p1 = mean(p1, 1);
p2 = mean(p2, 1);
p3 = mean(p3, 1);
p4 = mean(p4, 1);
p5 = mean(p5, 1);
p6 = mean(p6, 1);
p7 = mean(p7, 1);
p8 = mean(p8, 1);

prs = [p1; p2; p3; p4; p5; p6;p7;p8];
pstds = [p1s; p2s; p3s; p4s; p5s; p6s;p7s;p8s];
end

function [train, semitrain, test] = get_l_ul_us_inx(L, l, ul)
c =unique(L);
NC = length(c);
train = [];
semitrain = [];
test = [];
for ii = 1 : NC
    inx = find(L == c(ii));
    s = sampling(1 : length(inx), l);
    train = [train;inx(s)];
    inx = setdiff(inx, inx(s));
    s = sampling(1 : length(inx), ul);
    semitrain = [semitrain;inx(s)];
    inx = setdiff(inx, inx(s));
    test = [test;inx];
end
end
function [train, semitrain, test] = get_l_ul_us_inx2(L, lratio, ulratio)
c =unique(L);
NC = length(c);
train = [];
semitrain = [];
test = [];
for ii = 1 : NC
    inx = find(L == c(ii));
    l = floor(length(inx) * lratio);
    ul = floor(length(inx) * ulratio);
    s = sampling(1 : length(inx), l);
    train = [train;inx(s)];
    inx = setdiff(inx, inx(s));
    s = sampling(1 : length(inx), ul);
    semitrain = [semitrain;inx(s)];
    inx = setdiff(inx, inx(s));
    test = [test;inx];
end
end
function s = sampling(R, n)
% 选择抽样，R为记录集合，n为抽取的样本数
% 算法参考：D. E. Knuth, TAOCP, vol.2, pp142，稍有改动
  
% 编写函数时用的测试数据
if ~nargin
    R = 1 : 8;
    n = 4;
end
  
N = length(R);
t = 0;   % 处理过的记录总数
m = 0;   % 已选得的记录数
  
while 1
    U  = rand;
    if (N-t)*U < n-m
        m = m + 1;
        s(m) = R(t+1);
        % 若已抽取到足够的记录，则算法终止
        if m >= n, break, end
    end
    t = t + 1;
end
end

function [XTraink, LTraink] = find_near_k_class(XTestkk, XTrain, LTrain, k)
distance = pdist2(XTestkk, XTrain);
mean_diss = [];
ultrain = unique(LTrain);
for ii = 1 : length(ultrain)
    mean_diss = [mean_diss;sum(sum(distance(:, LTrain == ultrain(ii))))];
end
[sorted,inxs] = sort(mean_diss, 'ascend');
XTraink = [];
LTraink = [];
for ii = 1 : k
    XTrainii = XTrain(LTrain == ultrain(inxs(ii)), :);
    XTraink = [XTraink;XTrainii];
    LTraink = [LTraink;LTrain(LTrain == ultrain(inxs(ii)))];
end
end


% for ii = 1 : random_times
%     %try
%     l =5;
%     ul =3;
%     us = 7;
%     [train, semitrain, test] = get_l_ul_us_inx(L, l, ul);
%     %[train, semitrain, test] = get_l_ul_us_inx2(L, 0.5, 0.1);
%     XLTrain = X(train, :);
%     LTrain = L(train);
%     XULTrain = X(semitrain, :);
%     XTest = X(test, :);
%     LTest = L(test);
%  
%     [eigvector, W] = semi_dimension_reduction_main(XLTrain, LTrain, XULTrain, ALg_type);
%     XLTrain1 = XLTrain * eigvector1 * W(:, 1 : sel_dim);
%     XULTrain1 = XULTrain * eigvector1 * W(:, 1 : sel_dim);
%     XTest1 = XTest * eigvector1 * W(:, 1 : sel_dim);  
%     
%    for jj = 1 : length(sel_dims)
%         str = sprintf('%d_%d', ii, jj);
%         filename1 = strcat(file_name, str);
%         
%         sel_dim = sel_dims(jj);
%         ultest = unique(LTest);
%         classA = [];
%         LTestA = [];
%         for kk = 1 : length(ultest)
%             XTestkk = XTest(LTest == ultest(kk), :);
%             LTestkk = LTest(LTest == ultest(kk));
%             [XTraink, LTraink] = find_near_k_class(XTestkk, XLTrain, LTrain, length(ultest)* 1/2);
%             [eigvector, W] = semi_dimension_reduction_main(XTraink, LTraink, XULTrain, ALg_type);
%             if sel_dim > size(W, 2)
%                sel_dim = size(W, 2);
%             end
%             XLTrain1 = XTraink * eigvector * W(:, 1 : sel_dim);
%             XULTrain1 = XULTrain * eigvector * W(:, 1 : sel_dim);
%             XTest1 = XTestkk * eigvector * W(:, 1 : sel_dim);  
%             [ classes ] = my1nn( XLTrain1, LTraink, XTest1);
%             classA = [classA;classes];
%             LTestA = [LTestA;LTestkk];
%         end
%         p = length(LTestA(LTestA == classA))/length(classA);
%         p1(ii, jj) = p;
%         
%         sel_dim = sel_dims(jj);
%         if sel_dim > size(W1, 2)
%             sel_dim = size(W1, 2);
%         end
%         XLTrain1 = XLTrain * eigvector1 * W1(:, 1 : sel_dim);
%         XULTrain1 = XULTrain * eigvector1 * W1(:, 1 : sel_dim);
%         XTest1 = XTest * eigvector1 * W1(:, 1 : sel_dim);  
%         %[ classes ] = LMC1( XTest1', XLTrain1', LTrain, 3);
%         [ classes ] = my1nn( XLTrain1, LTrain, XTest1);
%         p = length(LTest(classes == LTest))/length(classes);
%         p2(ii, jj) = p;
%         
%         if sel_dims(jj) > size(W, 2)
%             break
%         end
%    end     
% end