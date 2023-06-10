function [ SF ] = get_graph_farL(  XL, L, K  )

X = XL';
[N, D] = size(X);

distance = pdist2(X, X);
[sorted,index] = sort(distance);
neighbor = zeros(N, N);
for ii = 1 :N
    neighbor(ii, index(ii, 1:K+1)) = 1;
end

%¸³ÖµÓÐlabelµÄ
SF = zeros(N, N);
for ii = 1 : N
    for jj = 1 : N
        if(L(ii) ~= L(jj)) && ii ~= jj && neighbor(ii, jj) == 0;
            SF(ii, jj) = 1;
        end
    end
     SF(ii, :) =  SF(ii, :) ./sum(SF(ii, :));
end

SF=SF-diag(diag(SF));
SF = (SF + SF')/2;

end

