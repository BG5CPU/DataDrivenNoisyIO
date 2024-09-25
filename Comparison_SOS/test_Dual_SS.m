% Implementation of 
% "Superstabilizing Control of Discrete-Time ARX Models under Error in Variables"
% Dual method

clear
clc
yalmip('clear')
% rng(1)

addpath('utils')
addpath('psatz')
addpath('Algorithm')


% define parameters
eps = [0.0001;0.0001]  ;          % noise bound    dy,du
d = 1;                        % degree of psatz
opts = sdpsettings('solver','mosek','verbose', 0);

% generate system
Adb1 = [cos(pi*0.5), sin(pi*0.5); -sin(pi*0.5), cos(pi*0.5)];
Adb2 = [cos(pi*0.4), sin(pi*0.4); -sin(pi*0.4), cos(pi*0.4)];
Adb3 = [cos(pi*0.3), sin(pi*0.3); -sin(pi*0.3), cos(pi*0.3)];
Adb4 = [cos(pi*0.2), sin(pi*0.2); -sin(pi*0.2), cos(pi*0.2)];
Adb5 = [cos(pi*0.1), sin(pi*0.1); -sin(pi*0.1), cos(pi*0.1)];
Ad = blkdiag(Adb2, Adb3, Adb4);
dimN = size(Ad,1);
dimM = 1;
dimP = 1;
Bd = ones(dimN,dimM);
Cd = ones(dimP,dimN);
Dd = 0;
[op_numb, op_denoa] = ss2tf(Ad,Bd,Cd,Dd);
Gz = tf(op_numb,op_denoa,0.1);


% transform to system in lambda
[Gl,al,bl] = sys_trans(Gz);     % al, bl (low to high order)
na_g = length(al);              % size of model a
nb_g = length(bl);              % size of model b
na_c = na_g+2;                       % size of controller a
nb_c = na_c-1;                       % size of controller b

n = [na_g;nb_g;na_c;nb_c];

T = 4*dimN-na_g; % # of samples

% generate trajectory
y = zeros(T+na_g,1);
y(1:na_g) = rand(na_g,1); 
u = (rand(T+nb_g-1,1)-0.5)*2*1;        

for i = 1:T
    y(i+na_g) = -al'*y(i+na_g-1:-1:i) + bl'*u(i+nb_g-1:-1:i);
end

y_noise = y + eps(1)*( 2*( rand(T+na_g,1)-0.5 ) );
u_noise = u + eps(2)*( 2*( rand(T+nb_g-1,1)-0.5 ) );

sim = struct('y_noise',y_noise,'u_noise',u_noise,'epsilon',eps);


tic

% solve Alternative
out = Alt_SS(sim, n, d, T);
sol = optimize(out.cons, out.obj, opts)

t = toc;

% extract solution
gamma = value(out.obj)
ac = value(out.ac);
bc = value(out.bc);

[Cz, ACL_z, Cl, ACL_l] = recover_sol(Gl, ac, bc)


disp(t);
disp(y_noise(1));

rmpath('utils')
rmpath('psatz')
rmpath('Algorithm')