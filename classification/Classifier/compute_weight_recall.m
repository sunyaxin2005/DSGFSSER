function [ weight_recall ] = compute_weight_recall( classes, LTest )
%COMPUTE_WEIGHT_RECALL Summary of this function goes here
%   Detailed explanation goes here
UL = unique(LTest);
weight_recall = 0;
for ii = 1 : length(UL)
    LTestii = LTest(LTest == UL(ii));
    classesii = classes(LTest == UL(ii));
    recallii = length(classesii(classesii == LTestii))/length(LTestii);
    weight_recall = weight_recall +  recallii * length(LTestii)/length(LTest);
end

end

