function [ unweight_recall ] = compute_unweight_recall( classes, LTest )
%COMPUTE_UNWEIGHT_ACCURY Summary of this function goes here
%   Detailed explanation goes here
%classes: recognition results
%LTest: class Label
UL = unique(LTest);
unweight_recall = 0;
for ii = 1 : length(UL)
    LTestii = LTest(LTest == UL(ii));
    classesii = classes(LTest == UL(ii));
    recallii = length(classesii(classesii == LTestii))/length(LTestii);
    %fprintf('%f ', recallii);
    unweight_recall = unweight_recall +  recallii;
end
unweight_recall = unweight_recall/length(UL);
end

