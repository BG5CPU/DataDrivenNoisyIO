
function [ZZ,index] = momentPowers(n1,n2,d)
% n: # of variables = n1+n2 (n1 is # of 0-1 variables, n2 is # of normal
% variables)
% d: degree of polynomial
% this file is part of SOSTOOLS, modified by NO
% generates the moments of at most degree d in R^n in lexicographical order
% and outputs an M(d) by n list containing powers
% here M(d) = nchoosek(n+d,n)

n = n1+n2;
ZZ = sparse(1,n);
for i = 1:n
    ss = size(ZZ,1);
    % ZZ = sprepmat(ZZ,d+1,1);
    ZZ = repmat(ZZ,d+1,1);
    for j = 0:d
        ZZ(ss*j+1:ss*j+ss,i) = j;
    end
%     idx = find(sum(ZZ,2) <= d);   % Throw away invalid monomials
    ZZ = ZZ(sum(ZZ,2)<=d,:);
end
ZZ = full(ZZ);

%% eliminate high order elements of n1 variables 
[mz,~] = size(ZZ);

elem = zeros(mz,1);
if n1 ~= 0
    for i = 1:mz
        temp = ZZ(i,:);
        temp(1:n1) = temp(1:n1)>=1;
%         rep = repmat(temp,mz,1);
%         equ = rep == ZZ;
        equ = bsxfun(@eq,ZZ,temp);
        num = find(sum(equ,2)==n);
        if i ~= num
            elem(i) = i;
        end
    end
end
elem(elem==0) = [];
ZZ(elem,:) = [];

%% range dictionary in the order of increasing degree
temp = cell(d+1,1);
for i = 0:d
    temp{i+1} = ZZ(sum(ZZ,2)==i,:);
end
ZZ = cell2mat(temp);

index = 1:size(ZZ,1);

end