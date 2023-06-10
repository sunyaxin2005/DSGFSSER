function G  = OLDA(Data, Class)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:                                                             %
%     Data:  data matrix (Each row is a data point)                  %
%     Class: class label (class 1, ..., k)                           %
% Output:                                                            %
%      G:    transformation matrix                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

k        = max(Class); % number of classes

%-------------------------------------------------------------------------
[m,n] = size(Data);
cc    = sum(Data)/m;

for i = 1:k 
	loc                  = find(Class==i);
	[num,o]              = size(loc);

	if num==1
		B(loc,1:n)   = 0;
		B(m+i,1:n)   = 0;
	else
		TMP          = sum(Data(loc,1:n))/num;
               
		B(i,1:n)     = sqrt(num)*(TMP - cc);
		B(loc+k,1:n) = Data(loc,1:n) - ones(num,1)*TMP;
	end;
end;

Hw = B(k+1:k+m, :)';

Hb = B(1:k, :)';

Size_low = rank(Hb);

Ht = (Data - ones(m,1)*cc)';

[U1,D1,V1] = svd(Ht,0);     % Ht = U1 D1 V1'
s = rank(Ht);
                                                                                                      
                                                                                                      
D1 = D1(1:s, 1:s);
U1 = U1(:,1:s);
                                                                                                      
d1 = diag(D1);
d1 = 1./d1;
D1 = diag(d1);

B = D1*U1'*Hb;

[P,D2, Q] = svd(B);

X = U1*D1*P;

S = X(:, 1:Size_low);

[Q,R] = qr(S,0);

G = Q;
