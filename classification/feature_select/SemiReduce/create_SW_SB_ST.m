function [ SW, SB, ST ] = create_SW_SB_ST( XL, L, XUL, W )
%CREATE_SW_SB Summary of this function goes here
%   Detailed explanation goes here
fea = [XL;XUL];
[nSmp,nFea] = size(fea);
D = full(sum(W,2));
W = -W;
for i=1:size(W,1)
    W(i,i) = W(i,i) + D(i);
end
sampleMean = mean(fea,1);
fea = (fea - repmat(sampleMean,nSmp,1));
D = W;
ST = fea'*D*fea;
ST = max(ST, ST');

classLabel = unique(L);
nClass = length(classLabel);
Hb = zeros(nClass,nFea);
for i = 1:nClass,
    index = find(gnd==classLabel(i));
    classMean = mean(XL(index,:),1);
    Hb (i,:) = sqrt(length(index))*classMean;
end
SB = Hb'*Hb;
SB = max(SB, SB');



end

