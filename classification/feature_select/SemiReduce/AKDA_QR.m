function [mean_accuracy, mean_time]  = AKDA_QR(dataX, dataY, k, rbf_var, F, reg)
% dataX:   data matrix stored in row-wise format
% dataY:   class label
% k:       number of nearest neighbors in KNN
% rbf_var: gaussian bandwidth  
% F:       number of training samples per class
% ref:     regularization parameter

mean_time = 0;

nnn = size(dataX,1);

av = mean(dataX);
sttd = std(dataX);

if 1==1
for i = 1:nnn
     dataX(i,:) = dataX(i,:)-av;
end
for i = 1:length(sttd)
    dataX(:,i) = dataX(:,i)/sttd(i);
end
end;

% dataX: data matrix stored in row-wise format
% dataY: class label  
% k:     number of nearest neighbors in KNN

classNum = max(dataY);  % number of classes

typeCollection = zeros(k,1);
weightscheme = 0; % set true

% ten pieces
for fold = 1: F
   pieceX{fold} = dataX(fold:F:size(dataX,1),:);
   pieceT{fold} = dataY(fold:F:size(dataY,1));
end

accuracy = zeros(F,1);

for round = 1: F

  trainingX = [];    % Training data
  trainingT = [];

  for id = 1: F
    if id ~= round
       trainingX = [trainingX; pieceX{id}];
       trainingT = [trainingT; pieceT{id}];
    end
  end
  testX = pieceX{round};    % Test data
  testT = pieceT{round};
 trainingNum = size(trainingX,1);
 testNum = size(testX,1);


 [V, R, centroids] = Kernel_QR_LDA(trainingX, classNum, rbf_var, 
trainingT, reg);

 trainingX_ = V'*inv(R')*RBFKernel(centroids,trainingX',rbf_var);
 trainingX  = trainingX_';
 
 testX_ = V'*inv(R')*RBFKernel(centroids,testX',rbf_var);
 testX  = testX_';
 
  sims = zeros(k,1);
  predict = zeros(testNum, 1);
  dists = zeros(trainingNum,1);

  % compute Mahalanobis kernal
%  Mahalanobis = pinv(trainingX' * trainingX);

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % begin to train
 for i= 1: testNum
   % compute the distances between ith test data and all training data
   for j = 1: trainingNum
%     sim = testX(i,:) * trainingX(j,:)' / (norm(testX(i,:)) * 
norm(trainingX(j,:)));
%     dists(j) = 1/ (1 + sim);
%      sim = testX(i,:) * Mahalanobis * trainingX(j,:)';
%     dists(j) = 1/ (1 + sim);
     dists(j) = norm(testX(i,:) - trainingX(j,:),2);
   end
  
   votes = zeros(classNum,1);

   % count the majority of classtype
   if ~weightscheme
    [orderedDists, I] = sort(dists);
    for l = 1: k
     typeCollection(l) = trainingT(I(l));
    end 
    for l = 1: k
       votes(typeCollection(l)) = votes(typeCollection(l)) + 1;
    end  
    [max_term id] = max(votes);
%    predict(i)= typeCollection(1);
    predict(i)= id;

  else % apply weight scheme
   
   [orderedDists, I] = sort(dists);
   for l = 1: k
     typeCollection(l) = trainingT(I(l));
     sims(l) = 1/ (1+ orderedDists(I(l)));
   end 
   probabilities = sims / sum(sims);  
    for l = 1: k
       votes(typeCollection(l)) = votes(typeCollection(l))  + sims(l);
    end  
    [max_term id] = max(votes);
    plot(votes); pause;
    predict(i)= id;

  end 

end

% count accuracy
diff = predict - testT;
idx = find(diff == 0);  % predict is correct
accuracy(round) = length(idx) / testNum;
%sss = length(idx) / testNum
accu = accuracy'
end

mean_accuracy = mean(accuracy)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [V, R, centroids] = Kernel_QR_LDA(trainingdata,classnum, rbf_var, 
classlabel, reg)

n = size(trainingdata,1);
k = classnum;
d = size(trainingdata,2);
av = mean(trainingdata)';

trainingdata = trainingdata';


centroids = [];

for s = 1:k
    loc = find(classlabel==s);
    
    x = mean(trainingdata(:,loc)')';
    
    centroids = [centroids, x];
end

CK = RBFKernel(centroids,centroids,rbf_var);
R = chol(CK);
Hb = centroids-av*ones(1,k);
CHb = RBFKernel(centroids,Hb,rbf_var);
Sb = inv(R')*CHb*CHb'*inv(R);


Ht = trainingdata-av*ones(1,n);
CHt = RBFKernel(centroids,Ht,rbf_var);
St = inv(R')*CHt*CHt'*inv(R);

[V,D] = eig(inv(St+ reg*eye(k))*Sb);

return




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function K = RBFKernel(X,Y,rbf_var)
xnum = size(X,2);
ynum = size(Y,2);
for i=1:xnum,
  for j=1:ynum,
    K(i,j) = exp(-norm(X(:,i)-Y(:,j))^2/rbf_var);
  end
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DP = DataProjection (V, R, C, rbf_var, testdata)
DP = V'*inv(R')*RBFKernel(C,testdata,rbf_var);
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


