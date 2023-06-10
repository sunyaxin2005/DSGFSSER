function [ output_args ] = MWMR( W, A )
%MWMR Summary of this function goes here
%   Detailed explanation goes here
[C,I] = max(W);
y = zeros(length(W), 1);
y(I) = 1;
inx = 1 : length(W);
for ii = 1 : k
    S1 = (y == 1);
    S3 = (y == 0);
    S2 = (y > 0 & y < 1);
    r = W/k - 2 * A *y/(k(k-1));
    
    F1 = S1 | S2;
    F2 = S2 | S3;
    rF1 = r(F1);
    rF2 = r(F2);
    inxF1 = inx(F1);
    inxF2 = inx(F2);
    
    [~,i] =max(rF1);
    [~,j] =max(rF2);
    i = inxF1(i);
    j = inxF1(j);
    
    CA = 2 * A(i, j) - A(i, i)-
    if r(i) > r(j)
        a = 
    
end
end

