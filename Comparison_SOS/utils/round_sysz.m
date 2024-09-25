function out = round_sysz(sys)
sys = minreal(sys);
[b,a] = tfdata(sys);
b = b{1};
a = a{1};
b = b/a(1);
a = a/a(1);
b(abs(b)<1e-5) = 0;
a(abs(a)<1e-5) = 0;
out = tf(b,a,0.1);
end