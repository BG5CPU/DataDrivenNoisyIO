function [vP, vK] = findControllerPetersenMIMO(dimN, dimP, dimM, Adz, Bdz, Cdz, A5, Bu, LL)

tune = 1e-15;

dimAll = dimN*dimP+dimN*dimM;

cvx_begin sdp
    variable P(dimAll, dimAll) semidefinite;
    variable Y(dimM, dimAll);

    blockS = [-P-LL*Cdz*LL',  A5*P+Bu*Y,  LL*Bdz;
               P*A5'+Y'*Bu',         -P,      -P;
               Bdz'*LL',             -P,    -Adz];
        
    minimize( -trace(P) );
    
    subject to   
        blockS <= 0;
        P >= tune*eye(dimAll);
cvx_end

vP = P; 
vK = Y/P;

end