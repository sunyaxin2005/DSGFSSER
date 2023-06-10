function [w,CPU_TIME]  = DRLDA(Data, Class)
% [w,CPU_TIME]  = DRLDA(Data, Class)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:                                                               %
%     Data:  data matrix (Each row is a data point)                    %
%     Class: class label (class 1, ..., c)                             %
% Output:                                                              %
%          w: transformation matrix                                    %
%  CPU_TIME : computational time                                       %
%                                                                      %
% Sharma, A and Paliwal, KK., A Deterministic Approach to Regularized  %
%    Linear Discriminant Analysis, Neurocomputing, 2014.                %
% For any comments or suggestions please email: alok.fj@gmail.com      %
%                                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
t0 = cputime;

%------------------------------ START -------------------------------%
c        = max(Class); % number of classes
[Nvec,dim] = size(Data);
mu    = sum(Data)/Nvec;

for i = 1:c
	loc     = find(Class==i);
	[nvec(i),col] = size(loc);
end

% STEP 1: Pre-processing stage
Ht = (Data - ones(Nvec,1)*mu)';
[Ut,Dt] = svd(Ht,0);
rt=rank(Dt);
Ut = Ut(:,1:rt);
m = Ut'*Data';
m=m';
mu = mean(m);
Hb=[]; Hw=[];
for i=1:c
    rg = 1+sum(nvec(1:i-1)):sum(nvec(1:i));
    x = m(rg,:);
    mui = mean(x);
    Hb = [Hb, sqrt(nvec(i))*(mui-mu)'];
    Hw = [Hw, (x - ones(nvec(i),1)*mui)'];
end
rb = rank(Hb);


% STEP 2: Find highest eigenvalue lambda_max
Sb = Hb*Hb';
Sw = Hw*Hw';
[Uw,Dw]=svd(Hw,0);
rw = rank(Dw);
Dw = Dw(1:rw,1:rw);
Uw = Uw(:,1:rw);
Swh = Uw*(Dw^-2)*Uw';
[wh,Dh] = eig(Swh*Sb);
dh = diag(Dh);
[Eval,inx] = sort(dh,'descend');
lambda_max = Eval(1); 

% STEP 3: Compute alpha
Q = (Sb/lambda_max - Sw);
[E,D] = eig(Q);
d = diag(D);
[alpha,inx] = sort(d,'descend');

% STEP 4: Compute orientation matrix W in rt-dimensional space
R = inv(Sw + alpha(1)*eye(rt))*Sb;
[E2,D2] = eig(R);
d2 = diag(D2);
[gamma,inx] = sort(d2,'descend');
E2 = E2(:,inx(1:rb));
w = E2;

% STEP 5: Compute orientation matrix in d-dimensional space
w = Ut*w;
%------------------------------  END  -------------------------------%

CPU_TIME = cputime - t0;
end
