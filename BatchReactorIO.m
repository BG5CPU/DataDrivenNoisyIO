clear; clc; close;

td    =  0.2; % sampling time
dimN  =  4;    % dimension of state
dimM  =  2;    % dimension of input
dimP  =  2;    % dimension of output
dimL  =  2;

%% parameter setting ======================================================
% system parameter setting 

Act = [ 1.38,   -0.2077, 6.715, -5.676;
       -0.5814, -4.29,   0,      0.675;
        1.067,   4.273, -6.654,  5.893;
        0.048,   4.273,  1.343, -2.104 ];
Bct = [ 0,      0;
        5.679,  0;
        1.136, -3.146;
        1.136,  0];
Cct = [1, 0, 1, -1;
       0, 1, 0,  0];
Dct = zeros(dimP,dimM);

Sct =  ss(Act, Bct, Cct, Dct);
Sdt =  c2d(Sct, td, 'zoh');
Ad  =  Sdt.A;
Bd  =  Sdt.B;
Cd  =  Sdt.C;
% Dd  =  Sdt.D;



%% calculate I/O matrix ===================================================
OB    =  [Cd' , (Cd*Ad)']';
Gama  =  [zeros(2), zeros(2); Cd*Bd, zeros(2)];
Rl    =  [Cd*Ad*Bd, Cd*Bd];
ioAB  =  [Cd*Ad^2/OB, Rl-Cd*Ad^2/OB*Gama];
F1    =  ioAB(:,1:2);
F2    =  ioAB(:,3:4);
G1    =  ioAB(:,5:6);
G2    =  ioAB(:,7:8);



%% generate data ==========================================================
Nsample =  5;
trajt   =  10;

eu = 20;    % bound of input
ed = 0.01;  % bound of output noise
ew = 0.01;  % bound of input noise

in0 = 2;
inT = 5;

U0 = []; X0 = [];
X1 = []; D0 = [];

for it = 1:trajt
    uu = (rand(dimM, Nsample)-0.5)*2* eu;
    xx = zeros(dimN, Nsample+1);
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


% plot(time, yd);
Bell  = [ zeros(dimP,dimM); zeros(dimP,dimM); zeros(dimM); eye(dimM) ];

A3ell = [ zeros(dimP), eye(dimP);
          zeros(dimP), zeros(dimP);
          zeros(dimM*dimL,dimP*dimL) ];

A4ell = [ zeros(dimP*dimL,dimM*dimL);
          zeros(dimM), eye(dimM); 
          zeros(dimM), zeros(dimM)];
Fell  = [ A3ell, A4ell ];
Lell  = [ zeros(dimP); eye(dimP); zeros(dimM,dimP); zeros(dimM,dimP)];

Aell  = [ [zeros(dimP), eye(dimP); F1, F2],... 
          [zeros(dimP,dimM), zeros(dimP,dimM); G1, G2];...
           zeros(dimM*dimL,dimP*dimL),...    
          [zeros(dimM), eye(dimM); zeros(dimM), zeros(dimM)] ];
Bdr   = [ zeros(dimP), zeros(dimP), zeros(dimP), zeros(dimP,dimM), zeros(dimP,dimM);
          -F1, -F2, eye(dimP), -G1, -G2;
          zeros(dimM*dimL,dimP*dimL+dimP+dimM*dimL) ];
% check the data, Tes should be zero
Tes =  Aell*X0 + Bell*U0 + Bdr*D0 - X1;


Ndat =  size(U0,2);
Dd   =  ones(dimP*dimL+dimP, Ndat)*ed;
Dw   =  ones(dimM*dimL, Ndat)*ew;
DAll =  [Dd;Dw];
De   =  max(svd(DAll*DAll')) * eye(dimP*dimL+dimP+dimM*dimL);

tm1 =  dimP;
tm2  = dimP*dimL+dimP+dimM*dimL;
De11 = De(1:tm1,     1:tm1); De12 = De(1:tm1,     tm1+1:tm2); 
De21 = De(tm1+1:tm2, 1:tm1); De22 = De(tm1+1:tm2, tm1+1:tm2); 

Adz  = X0*X0'-De22;
Bdz  = -Lell'*X1*X0' + De12;
Cdz  = (Lell'*X1)*(Lell'*X1)' - De11;

eig_X0  = eig(X0*X0');
eigDe22 = eig(De22);
eigAdz  = eig(Adz);




%% find the controller ====================================================
[vP, vK] = findControllerPetersenMIMO(dimL, dimP, dimM, Adz, Bdz, Cdz, Fell, Bell, Lell);
vK
eig_vK = eig(Aell+Bell*vK)

close_FLB = [Ad,      Bd*vK;
             Lell*Cd, Fell+Bell*vK];

eig_close_FLB = eig(close_FLB)








%% state space realization of the controller===============================
K1y = vK(:,1:2); K2y = vK(:,3:4); 
K1u = vK(:,5:6); K2u = vK(:,7:8);

Ak = [zeros(2), K1u;
      eye(2),   K2u ];
Bk = [K1y; K2y];
Ck = [zeros(2), eye(2)];



%% state space realization of the closed-loop system=======================
Acl = [Ad,    Bd*Ck;
       Bk*Cd, Ak    ];

eig_Acl = eig(Acl);


 

















