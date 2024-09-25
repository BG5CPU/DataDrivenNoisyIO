function out = Alt_SS(sim, n, d, T)
%% Alternatives method, Algorithm 3
%  sim:     sampled trajectories
%    n:     size of coefficients
%    d:     degree of the psatz
%    T:     # of samples for design

% extract data and size
na_g = n(1);
nb_g = n(2);
na_c = n(3);
nb_c = n(4);
y_noise = sim.y_noise(1:T+na_g);        % noisy output
u_noise = sim.u_noise(1:T+nb_g-1);      % noisy input
eps = sim.epsilon;                      % noise bound epsilon
Cons = [];                              % save constraints

% define variables
ag = sdpvar(na_g,1);
bg = sdpvar(nb_g,1);
ac = sdpvar(na_c,1);
bc = sdpvar(nb_c,1);
vars = [ag;bg];

% define constraints
g = [eps(1)*ones(2*(T+na_g),1);eps(2)*ones(2*(T+nb_g-1),1)];   
for i = 1:T                  
    h(i,1) = y_noise(i+na_g) + ag'*y_noise(i+na_g-1:-1:i) - bg'*u_noise(i+nb_g-1:-1:i);
end
cons_data = struct('ineq', g, 'eq', h);    % constraints from data

% define m, acl
da = na_g+na_c; 
db = nb_g+nb_c;
for i = 1:da        % eq.(28c)
    m(i,1) = polynomial([ag;bg],2*d);
end
AA = compute_coeff([1;ag],[1;ac]);
BB = compute_coeff([0;bg],[0;bc]);
AA(1) = [];
BB(1) = [];
if length(AA) > length(BB)
    BB = [BB; zeros(da-db,1)];
end
acl = AA + BB;

% define non-negative polynomial eq.(28de)
gamma = sdpvar;
p1 = m-acl;
p2 = m+acl;
p3 = gamma-sum(m);
P = [p1(:);p2(:);p3(:)];    % P contains m-acl >= 0, m+acl >= 0, gamma-tol-sum(M,2) >= 0
l = length(P);

% get psatz constraints
Gram = cell(l,1);                       % save Gram matrix
Coef = cell(l,1);                       % save coef of mu
for i = 1:l
    [p_psatz, cons_psatz, Gram{i,1}, Coef{i,1}] = Dual_SS_psatz(P(i), cons_data, d, vars, n, T);
    Cons = [Cons; cons_psatz];          % constraints from psatz
end
out.poly = p_psatz;
out.cons = Cons;
out.Gram = Gram;
out.Coef = Coef;
out.obj = gamma;
out.ac = ac;
out.bc = bc;
end
