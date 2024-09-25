function out = compute_coeff(a,b)
% compute coefficients of a*b
% by sum over anti-diagonal elements
temp = a*b';
[m,n] = size(temp);
ind = hankel(1:m,m:(n-1)+m);
out = accumarray(ind(:),temp(:));
end



