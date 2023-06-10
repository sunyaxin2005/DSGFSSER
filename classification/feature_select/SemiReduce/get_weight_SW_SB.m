function [ SB, SW, ST ] = get_weight_SW_SB(  XL, L )
%GET_SWSB Summary of this function goes here
%   Detailed explanation goes here
%XL:有标签训练样本，维度×样本个数
UL = unique(L);
K = 10000;
for jj = 1 : length(UL)
   if ( K > length(L(L == jj)))
       K = length(L(L == jj));
   end
end
X = XL';
[N, D] = size(X);
distance = pdist2(X, X);
[sorted,index] = sort(distance);
sita = sorted(K+1, :);
%sita = sita' *  sita;
t = mean(mean(distance));
A = exp(-(distance.^2).^2./(t));
% diagA = diag(ones(size(distance, 2), 1));
% A = A - diagA;
% neighborhood = index(2:(1+K),:);
%A = ones(N,N);
sita = sita' *  sita;
distance  = distance - sqrt(sita);
distance(distance < 0) = 0;
%sita1 = sorted(2 * K+1, :);
% sita1 = sita1' *  sita1;
% sita1 = sita1 - sita;
t = mean(mean(distance));
A1 = distance.^2 * exp(-(distance.^2).^2./(t));
SW = zeros(N, N);
UL = unique(L);
nClass = length(UL);
%赋值有label的
% for ii = 1 : nClass
%     SW(L == UL(ii), L == UL(ii)) = A(L == UL(ii), L == UL(ii))/length(L(L == UL(ii)));
% end% %赋值没有label的
for ii = 1 : N
    for jj = 1 : N
        if(L(ii) == L(jj)) && ii ~= jj
            SW(ii, jj) = 1;
        end
    end
    SW(ii, :) = SW(ii, :)./sum(SW(ii, :));
end
SW = SW + SW'- SW' * SW;
%赋值没有label的
% for ii = 1 : N
%     SW(neighborhood(:, ii), ii) = SW(neighborhood(:, ii), ii) + A(neighborhood(:, ii), ii)/N;
%     SW(ii, neighborhood(:, ii)) = SW(ii, neighborhood(:, ii)) + A(ii, neighborhood(:, ii))/N;
% end
SW=SW-diag(diag(SW));
%赋值有label的
SB = zeros(N, N);
% for ii = 1 : nClass
%     SB(L == UL(ii), L == UL(ii)) = -A(L == UL(ii), L == UL(ii))*(1/length(L(L == UL(ii))) - 1/N);
%     SB(L == UL(ii), L ~= UL(ii)) = A(L == UL(ii), L ~= UL(ii))/N;
%     SB(L ~= UL(ii), L == UL(ii)) = A(L ~= UL(ii), L == UL(ii))/N;
% end% %赋值没有label的
for ii = 1 : N
    for jj = 1 : N
        if(L(ii) ~= L(jj)) && ii ~= jj
            SB(ii, jj) = A1(ii, jj);
        end
    end
     SB(ii, :) =  SB(ii, :) ./sum(SB(ii, :));
end
SB = SB + SB' - SB'*SB;
SB=SB-diag(diag(SB));
ST = A;
end



