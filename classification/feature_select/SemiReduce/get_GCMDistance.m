function [ Dis ] = get_GCMDistance( XTest, XTrain, K )
%GET_GCMDISTANCE Summary of this function goes here
%   Detailed explanation goes here

distance = pdist2(XTrain, XTrain);
[sorted,index] = sort(distance, 2);
sita1 = sorted(:, floor(2 * K)+1);
sita = mean(mean(sita1));
FDX = exp((-distance.^2)./(sita' * sita));
FDX = sum(FDX, 2);
PX = FDX./sum(FDX);
% 
IXTrain = -log(PX);
%UL = unique(LTrain);
% mean_C_dis = zeros(length(UL), 1);
% mean_IXTrain = zeros(length(UL), 1);
% IXTrain1 = IXTrain * IXTrain';
% for ii = 1 : length(UL)
%     distance2 = distance.^2;
%     mean_C_dis(ii) = sum(sum(distance2(LTrain == UL(ii), LTrain ~= UL(ii))));
%     mean_C_dis(ii) = mean_C_dis(ii) + sum(sum(distance2(LTrain ~= UL(ii), LTrain == UL(ii))));
%     mean_IXTrain(ii) = mean_IXTrain(ii) + sum(sum(IXTrain1(LTrain == UL(ii), LTrain ~= UL(ii))));
%     mean_IXTrain(ii) =  mean_IXTrain(ii) + sum(sum(IXTrain1(LTrain ~= UL(ii), LTrain == UL(ii))));    
% end
% 
distance = pdist2(XTest, XTrain);
[sorted,index] = sort(distance, 2);
sita2 = sorted(:, floor(2 * K)+1);
sita = mean(mean(sita2));
FDX = exp(-(distance.^2)./(sita' * sita));
FDX = sum(FDX, 2);
PX = FDX./sum(FDX);
IXTest = -log(PX);
IX = IXTest * IXTrain';
%fprintf('dim:%d\n', size(XTest, 2));
%fprintf('IX  :%5.3f %5.3f %5.3f  %5.3f\n', min(IX(:)), max(IX(:)), mean(IX(:)), std(IX(:)));
sita = sita2 * sita1';
%fprintf('sita:%5.3f %5.3f %5.3f  %5.3f  %5.3f\n', min(sita(:)), max(sita(:)), mean(sita(:)), std(sita(:)), std(IX(:))/ std(sita(:)));
%fprintf('sita:%5.3f %5.3f %5.3f  %5.3f  %5.3f\n', min(distance(:)), max(distance(:)), mean(distance(:)), std(distance(:)), std(sita(:))/ std(distance(:)));

if std(IX(:)) ~= 0
    IX = (IX  - mean(mean(IX))) * 1.4 * std(sita(:))/(std(IX(:))) +  mean(mean(sita));
end
IX(IX < min(min(sita))) = min(min(sita));
PPP = -(distance.^2)./IX;
try
Dis = exp(PPP);
catch
    c= 1;
end

% Dis = distance./sqrt(IX);
% % for ii = 1 :size(Dis, 1)
% %     Dis(ii, :) = (Dis(ii, :) - min(Dis(ii, :)))./(max(Dis(ii, :))-min(Dis(ii, :)));
% % end
% [sorted,index] = sort(distance, 2);
% sita1 = sorted(:, 2 * K+1);
% sita = mean(mean(sita1));
% Dis = exp(-(Dis.^2)./(sita * sita));
% % 
Dis(isinf(Dis)) = 0;
end


