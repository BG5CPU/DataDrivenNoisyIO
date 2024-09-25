function [p_psatz, cons_psatz, Gram, Coef] = Dual_SS_psatz(p, C, d, vars, n, T)
%% Enforce Putinar Psatz using Algorithm 2
%    p:   nonnegative polynomial
%    C:   eq and ineq constraints
%    d:   degree of p
% vars:   variables a,b
%    n:   size of coefficients
%    T:   # of samples

na = n(1);
nb = n(2);
sy = T+na;                  % size of dy
su = T+nb-1;                % size of du
n_var = length(vars);       % # of variables
d_g = degree(C.ineq);       % degree of inequality
d_h = degree(C.eq);         % degree of equality
Gram = cell(1+2*(sy+su), 1);      % gram coefficients of sigma0, zeta, psi
Coef = cell(T, 1);          % coefficients of mu
p_psatz = 0;                % Putinar form of p

zeta_p = zeros(sy,1,'like',sdpvar);       % zeta corresponding to dy
zeta_n = zeros(sy,1,'like',sdpvar);
psi_p = zeros(su,1,'like',sdpvar);        % psi corresponding to du
psi_n = zeros(su,1,'like',sdpvar);
mu = zeros(T,1,'like',sdpvar);

% generate -Q
[pow0, ~] = momentPowers(0, n_var, d);  
vect_0 = recovermonoms(pow0, vars);       % basis v(x)
l_0 = length(vect_0);
Gram{1} = sdpvar(l_0);
Q = -vect_0'*Gram{1}*vect_0;
p_psatz = p_psatz - Q;
cons_psatz = (Gram{1} >= 0):'Gram_0';

% generate zeta, psi eq.(27ab)
if mod(d_g,2) == 0                        % decide degree of poly
    [pow_i, ~] = momentPowers(0, n_var, floor((2*d-d_g)/2));
else
    [pow_i, ~] = momentPowers(0, n_var, floor((2*d-d_g-1)/2));
end
vect_i = recovermonoms(pow_i, vars);
l_i = length(vect_i);

for i = 1:sy   % zeta+
    Gram{i+1} = sdpvar(l_i);
    zeta_p(i) = vect_i'*Gram{i+1}*vect_i;
    p_psatz = p_psatz + zeta_p(i)*C.ineq(i);
    cons_psatz = [cons_psatz; (Gram{i+1} >= 0):'Gram_zeta_p'];
end

for i = 1:sy    % zeta-
    Gram{i+1+sy} = sdpvar(l_i);
    zeta_n(i) = vect_i'*Gram{i+1+sy}*vect_i;
    p_psatz = p_psatz + zeta_n(i)*C.ineq(i+sy);
    cons_psatz = [cons_psatz; (Gram{i+1+sy} >= 0):'Gram_zeta_n'];
end

for i = 1:su   % psi+
    Gram{i+1+2*sy} = sdpvar(l_i);
    psi_p(i) = vect_i'*Gram{i+1+2*sy}*vect_i;
    p_psatz = p_psatz + psi_p(i)*C.ineq(i+2*sy);
    cons_psatz = [cons_psatz; (Gram{i+1+2*sy} >= 0):'Gram_psi_p'];
end

for i = 1:su     % psi-
    Gram{i+1+2*sy+su} = sdpvar(l_i);
    psi_n(i) = vect_i'*Gram{i+1+2*sy+su}*vect_i;
    p_psatz = p_psatz + psi_n(i)*C.ineq(i+2*sy+su);
    cons_psatz = [cons_psatz; (Gram{i+1+2*sy+su} >= 0):'Gram_psi_n'];
end

% generate mu eq.(27c)
[pow_mu, ~] = momentPowers(0, n_var, 2*d-d_h);
vect_mu = recovermonoms(pow_mu, vars);
l_mu = length(vect_mu);
for i = 1:T     
    Coef{i} = sdpvar(l_mu,1);
    mu(i) = vect_mu'*Coef{i};
    p_psatz = p_psatz + mu(i)*C.eq(i);
end

% eq.(27d-f)
cons_psatz = [cons_psatz; (coefficients(p - p_psatz, vars) == 0):'eq'];
cons_psatz = [cons_psatz; (coefficients(zeta_p - zeta_n - compute_corr(mu, [1;vars(1:na)]),vars) == 0):'eq_zeta'];
cons_psatz = [cons_psatz; (coefficients(psi_p - psi_n - compute_corr(mu, vars(na+1:end)),vars) == 0):'eq_psi'];
end


