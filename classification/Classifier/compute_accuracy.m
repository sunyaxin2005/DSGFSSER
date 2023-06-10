function [ accuracy ] = compute_accuracy( classes, LTest )
%COMPUTE_ACCURACY Summary of this function goes here
%   Detailed explanation goes here
accuracy = length(LTest(classes == LTest))/length(LTest);

end

