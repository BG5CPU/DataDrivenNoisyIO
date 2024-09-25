close all; clear; clc;


%% set state A B C D matrix ===============================================
% original system
Adb1 = [cos(pi*0.5), sin(pi*0.5); -sin(pi*0.5), cos(pi*0.5)];
Adb2 = [cos(pi*0.4), sin(pi*0.4); -sin(pi*0.4), cos(pi*0.4)];
Adb3 = [cos(pi*0.3), sin(pi*0.3); -sin(pi*0.3), cos(pi*0.3)];
Adb4 = [cos(pi*0.2), sin(pi*0.2); -sin(pi*0.2), cos(pi*0.2)];
Adb5 = [cos(pi*0.1), sin(pi*0.1); -sin(pi*0.1), cos(pi*0.1)];
Ad = blkdiag(Adb2, Adb3, 1);
dimN = size(Ad,1);
dimM = 1;
dimP = 1;
dimL = dimN;
Bd = ones(dimN,dimM);
Cd = ones(dimP,dimN);
Dd = 0;
[op_numb, op_denoa] = ss2tf(Ad,Bd,Cd,Dd);



%% generate data ==========================================================
Nsample =  4*dimN;
trajt   =  1;

eu = 1;    % bound of input
ed = 0.0001;  % bound of output noise
ew = 0.0001;  % bound of input noise

in0 = 1;
inT = Nsample;

U0 = []; X0 = [];
X1 = []; D0 = [];

for it = 1:trajt
    uu = (rand(dimM, Nsample)-0.5)*2* eu;
    xx = zeros(dimN, Nsample+1); xx(:,1) = randn(dimN,1)*0.5;
    yy = zeros(dimP, Nsample);
    dd = (rand(dimP, Nsample)-0.5)*2* ed;
    dw = (rand(dimM, Nsample)-0.5)*2* ew;    
    for ii = 1:Nsample
        yy(:, ii)   =  Cd*xx(:, ii);
        xx(:, ii+1) =  Ad*xx(:, ii) + Bd*uu(:, ii);
    end
    yd = yy+dd;
    uw = uu+dw;
    % creat data matrices
    U0_ = creatHankel(uw(:,in0+dimL:inT),1);
    X0_ = [creatHankel(yd(:,in0:inT-1),dimL) ; creatHankel(uw(:,in0:inT-1),dimL)];
    X1_ = [creatHankel(yd(:,in0+1:inT),dimL) ; creatHankel(uw(:,in0+1:inT),dimL)];
    D0_ = [creatHankel(dd(:,in0:inT),dimL+1) ; creatHankel(dw(:,in0:inT-1),dimL)];    
    U0  = [U0, U0_];
    X0  = [X0, X0_];
    X1  = [X1, X1_];
    D0  = [D0, D0_];
end


% noise energy bound
Theta = func_build_eTheta(ed, ew, 0, dimP, dimM, dimL, size(X0,2));
tm1 = dimP;
tm2 = dimP*dimL+dimP+dimM*dimL;
Theta11 = Theta(1:tm1,     1:tm1); Theta12 = Theta(1:tm1,     tm1+1:tm2); 
Theta21 = Theta(tm1+1:tm2, 1:tm1); Theta22 = Theta(tm1+1:tm2, tm1+1:tm2); 


Lell = func_build_Lell(dimP, dimM, dimL);
Bell = func_build_Bell(dimP, dimM, dimL);
Fell = func_build_Fell(dimP, dimM, dimL);


AdZ  = X0*X0'-Theta22;
BdZ  = -Lell'*X1*X0' + Theta12;
CdZ  = (Lell'*X1)*(Lell'*X1)' - Theta11;

eig_AdZ = eig(AdZ);
eig_AdZ_min = min(eig_AdZ);
disp("eig_AdZ_min is " + eig_AdZ_min);






%% find the controller ====================================================
tic

[vP, vK] = findControllerPetersenMIMO(dimL, dimP, dimM, AdZ, BdZ, CdZ, Fell, Bell, Lell);

t = toc;

disp(vK);
eig_vK = eig( [Ad, Bd*vK; Lell*Cd, Fell+Bell*vK] );
abs(eig_vK)
 
disp(t);











%% useful function ========================================================
function Fell = func_build_Fell(p, m, ell)
    blk1 = [ zeros( p*(ell-1), p ), eye( p*(ell-1) ); zeros( p, p*ell ) ];
    blk2 = zeros( p*ell, m*ell );
    blk3 = zeros( m*ell, p*ell );
    blk4 = [ zeros( m*(ell-1), m ), eye( m*(ell-1) ); zeros( m, m*ell ) ];
    Fell = [blk1, blk2; blk3, blk4];
end


function Lell = func_build_Lell(p, m, ell)
    Lell = zeros( p*ell+m*ell, p ); 
    Lell( p*(ell-1)+1 : p*ell, : ) = eye(p);
end


function Bell = func_build_Bell(p, m, ell)
    Bell = zeros( p*ell+m*ell, m ); 
    Bell( p*ell+m*(ell-1)+1 : p*ell+m*ell, : ) = eye(m);
end


function eTheta = func_build_eTheta(eY, eU, eYa, p, m, ell, clo)
    Dd  = ones(p*ell+p, clo)*eY;
    Dw  = ones(m*ell, clo)*eU;
    Da  = ones(p*ell+p, clo)*eYa;
    DAll = [Dd+Da;Dw];
    eTheta = max(svd(DAll*DAll')) * eye(p*ell+p+m*ell);
end




