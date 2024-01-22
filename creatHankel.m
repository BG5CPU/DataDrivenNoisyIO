function [mHankel] = creatHankel(seq,depL)

[m,n] = size(seq);
mHankel = zeros( m*depL , n-depL+1 );
for il = 1:depL
    st = (il-1)*m+1; ed = il*m;
    mHankel( st:ed , : ) = seq( : , il:n-depL+il );
end

end