function out = compute_corr(a,b)
% compute cross-correlation of a \star b
% by sum over off-diagonal elements
temp = a*b';
[m,n] = size(temp);
ind = hankel(1:m,m:(n-1)+m);
ind = flip(ind,2);
out = accumarray(ind(:),temp(:));
end

