function [ percision ] = computer_percision( classes, LTest )
%COMPUTER_PERCISION Summary of this function goes here
%   Detailed explanation goes here
%classes:n¡Á1,recognition resultes
%LTest:n¡Á1,the labels of testing points 
UL = unique(classes);
percision = 0;
for ii = 1 : length(UL)
    classesii = classes(classes == UL(ii));%get the recognition resultes of recogniting to UL(ii) class
    LTestii = LTest(classes == UL(ii));%get the corrsponding class labels of recogniting to UL(ii) class
    percisionii = length(classesii(classesii == LTestii))/length(classesii);
    percision = percision +  percisionii;
    %fprintf('%f ', percisionii)
end
percision = percision/length(UL);
end

