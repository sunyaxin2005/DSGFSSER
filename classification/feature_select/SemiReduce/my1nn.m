function [ L ] = my1nn( XLTrain, LTrain, XTest)
%MY1NN Summary of this function goes here
%   Detailed explanation goes here
distance = pdist2(XTest, XLTrain, 'euclidean');
[C, I] = min(distance, [], 2);
L = LTrain(I);

end

