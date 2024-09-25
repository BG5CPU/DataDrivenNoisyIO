function [sys2,a2,b2] = sys_trans(sys1)
% transfer between system in z and lambda, where lambda = 1/z
% only works if sys1 is strictly proper
[b1,a1] = tfdata(sys1);
b1 = b1{1};
a1 = a1{1};
if a1(1) ~= 0
a2 = flip(a1)/a1(1);
b2 = flip(b1)/a1(1);
else 
    a2 = flip(a1);
    b2 = flip(b1);
end
sys2 = tf(b2,a2,0.1);
a2(end) = [];
a2 = flip(a2)';
b2(end) = [];
b2 = flip(b2)';
while length(b2) > 1 & b2(end) == 0
    b2(end) = [];
end  
end