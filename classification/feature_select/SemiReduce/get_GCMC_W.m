function [ M ] = get_GCMC_W(  X, K )
%GET_GCMC_W Summary of this function goes here
%   Detailed explanation goes here
%distance = pdist2(X, X);
[W] = GCMDistance(X, X, K);
%W = distance .* W;
for ii=1:size(W, 2)
   W(:,ii) = W(:,ii)/sum(W(:,ii));                  % enforce sum(w)=1
end;
M = (W+W')/2;

end

function [Wg, Wl] = GetWlWg( XL, L, X )
UL = unique(L);
K = 10000;
for jj = 1 : length(UL)
   if ( K > length(L(L == jj)))
       K = length(L(L == jj));
   end
end
K = 2 * K;
distance = pdist2(X, X);
[sorted,index] = sort(distance);
sita = sorted(K+1, :);
sita = sita' *  sita;
t = mean(mean(sita));
Wg = distance.^2 * exp(-distance.^2./(2 * t.^2));
distance = pdist2(XL, XL);
[sorted,index] = sort(distance);
sita = sorted(K+1, :);
sita = sita' *  sita;
t = mean(mean(sita));
Wl = exp(-distance.^2./(2 * t.^2));
for ii = 1 : size(Wl, 1)
    for jj = 1 : size(Wl, 2)
        if(L(ii) ~= L(jj))
            Wl(ii, jj) = 0;
        end
    end
end
end