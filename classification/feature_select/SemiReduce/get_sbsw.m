function [ sb, sw ] = get_sbsw( XL, L )
[N, D] = size(XL);
UL = unique(L);
NC = length(UL);
miu = zeros(NC, D);
miuA = mean(XL, 1);
sb = zeros(D, D);
sw = zeros(D, D);
for ii = 1 : NC
    Xii = XL(L== UL(ii), :);
    miui = mean(Xii,1);
    miu(ii, :) = miui;
    sb = sb + size(Xii, 1) * (miui-miuA)' * (miui-miuA);
    for jj = 1 : size(Xii, 1)
        sw = sw + (Xii(jj, :) - miui)' * (Xii(jj, :) - miui);
    end
end


end

