function [ FDRs ] = FDR( X, L )
%FDR Summary of this function goes here
%   Detailed explanation goes here
UL = unique(L);
mean_ls = [];
std_ls = [];
for inx = 1 : length(UL)
    mean_ls = [mean_ls;mean(X(L == UL(inx), :))];
    std_ls = [std_ls;std(X(L == UL(inx), :)).^2];
end
FDRs = 0;
for inx1 = 1 : length(UL)
    for inx2 = inx1:length(UL)
        if (std_ls(inx1, :) + std_ls(inx2, :)) == 0
            FDRs = 0;
        else
            FDRs = FDRs + (mean_ls(inx1, :) - mean_ls(inx2, :)).^2./ (std_ls(inx1, :) + std_ls(inx2, :));
        end
        
    end
end
FDRs = FDRs * 2/((length(UL) - 1) * length(UL));
end

