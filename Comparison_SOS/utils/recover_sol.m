function [Cz, ACL_z, Cl, ACL_l] = recover_sol(Gl, ac, bc)
na_c = length(ac);
nb_c = length(bc);
ac = flip([1;ac])';     % add constant 1 and flip to high to low order
bc = flip([0;bc;zeros(na_c-nb_c,1)])';  % add and pad 0, flip

Cl = tf(bc,ac,0.1);     % controller in lambda
Cl = round_sysl(Cl);     % round off to drop small coefficients
Cz = sys_trans(Cl);     % controller in z
Cz = round_sysz(Cz);
ACL_l = Gl/(1+Gl*Cl);       % closed-loop in lambda
ACL_l = round_sysl(ACL_l);
ACL_z = sys_trans(ACL_l);   % closed-loop in z
ACL_z = round_sysz(ACL_z);
end
