  close all; clear; clc;

%% set state A B C D matrix ===============================================
% original system
Aini = [0, 1, 0;
       -1, 0, 0;
        0, 0, 1];
Bini = [1, 0; 0, 1; 1, 0];
Cini = [1, 0, 1; 0, 1, 1];
% artificial system
Aa   = 0;
Ba   = [1, 1];
Ca   = [1; 1];
% dimN=3 is not an integral multiple of dimP=2
% so we augment the stat
Aaug = blkdiag(Aini, 0.0);
Baug = [Bini; Ba];
Caug = [Cini, Ca];
r_ctrb = rank( ctrb(Aaug,Baug) );
r_obsv = rank( obsv(Aaug,Caug) );
OB = [Caug; (Caug*Aaug)];
rankOB = rank(OB);
% the augmented system is controllable and observable

dimP  = 2;  % dimension of input y
dimM  = 2;  % dimension of input u
dimN  = 3;  % dimension of original state 
dimNa = 1;  % dimension of artificial state 
dimAg = 4;  % dimension of augmented state 
dimL  = 2;


%% calculate I/O matrix ===================================================
OB    =  [Caug; (Caug*Aaug)];
Gama  =  [zeros(2), zeros(2); Caug*Baug, zeros(2)];
Rl    =  [Caug*Aaug*Baug, Caug*Baug];
ioAB  =  [Caug*Aaug^2/OB, Rl-Caug*Aaug^2/OB*Gama];
F1    =  ioAB(:,1:2);
F2    =  ioAB(:,3:4);
G1    =  ioAB(:,5:6);
G2    =  ioAB(:,7:8);


%% compare the output =====================================================
% % sysNM = drss(dimN,dimP,dimM); % generate random discrete test model
% 
% Nd = 200;             % data length
% 
% us = randn(dimM,Nd);  % pre-define the input sequence
% xs = zeros(dimN,Nd+1);
% ys = zeros(dimP,Nd);
% for ik = 1:Nd
%     xs(:,ik+1) = Aaug*xs(:,ik) + Baug*us(:,ik);
%     ys(:,ik)   = Caug*xs(:,ik);
% end
% figure(1); plot(ys(1,:), 'r'); hold on;
% 
% 
% yy = zeros(dimP,Nd+3);
% uu = [zeros(dimM,2), us];
% for ik = 3:Nd+3
%     yy(:,ik) = F2*yy(:,ik-1) +F1*yy(:,ik-2) +G2*uu(:,ik-1) +G1*uu(:,ik-2);
% end
% plot(yy(1,3:202), 'b:');




%% generate data ==========================================================
Ndata   =  30;
nStar   =  5;
Nsample =  40;

eu = 2;      % bound of input
ed = 0.01;   % bound of output noise
ew = 0.01;   % bound of input noise

% bound of output noise for artificial system
ea = 1* svd(Ca)* svd(Ba)* ew; 


uu = (rand(dimM, Nsample)-0.5)*2* eu;  % input
dw = (rand(dimM, Nsample)-0.5)*2* ew;  % input noise
dd = (rand(dimP, Nsample)-0.5)*2* ed;  % output noise


% data for the augmented system
xx = zeros(dimN, Nsample+1);
yy = zeros(dimP, Nsample);
for ii = 1:Nsample
    yy(:, ii)   =  Cini*xx(:, ii);
    xx(:, ii+1) =  Aini*xx(:, ii) + Bini*uu(:, ii);
end
xa = zeros(dimAg-dimN, Nsample+1);
ya = zeros(dimP, Nsample);
for ii = 1:Nsample
    ya(:, ii)   =  Ca*xa(:, ii);
    xa(:, ii+1) =  Aa*xa(:, ii) + Ba*(uu(:, ii)+dw(:,ii));
end
yda = (yy+dd) + ya;
uwa = uu+dw;


in0 = 2;
inT = in0+Ndata+1;
% creat data matrices
U0 = creatHankel(uwa(:,in0+dimL:inT),1);
X0 = [creatHankel(yda(:,in0:inT-1),dimL) ; creatHankel(uwa(:,in0:inT-1),dimL)];
X1 = [creatHankel(yda(:,in0+1:inT),dimL) ; creatHankel(uwa(:,in0+1:inT),dimL)];
D0 = [creatHankel(dd(:,in0:inT),dimL+1) ; creatHankel(dw(:,in0:inT-1),dimL)];

eig_X0 = eig(X0*X0');

Bell  = [ zeros(dimP,dimM); zeros(dimP,dimM); zeros(dimM); eye(dimM) ];

A3ell = [ [zeros(dimP), eye(dimP); zeros(dimP), zeros(dimP)]; 
          zeros(dimM*dimL,dimP*dimL) ];

A4ell = [ zeros(dimP*dimL,dimM*dimL);
          zeros(dimM), eye(dimM); 
          zeros(dimM), zeros(dimM) ];

Fell  = [ A3ell, A4ell ];

Lell  = [ zeros(dimP); eye(dimP); zeros(dimM,dimP); zeros(dimM,dimP) ];

Aell  = [ [zeros(dimP), eye(dimP); F1, F2],... 
          [zeros(dimP,dimM), zeros(dimP,dimM); G1, G2];...
           zeros(dimM*dimL,dimP*dimL),...    
          [zeros(dimM), eye(dimM); zeros(dimM), zeros(dimM)] ];

% Bdr   = [ zeros(dimP), zeros(dimP), zeros(dimP), zeros(dimP,dimM), zeros(dimP,dimM);
%           -F1, -F2, eye(dimP), -G1, -G2;
%           zeros(dimM*dimL,dimP*dimL+dimP+dimM*dimL) ];
%
% Tes =  Aell*X0 + Bell*U0 + Bdr*D0 - X1;

Ndat =  size(U0,2);
Dd   =  ones(dimP*dimL+dimP, Ndat)*ed;
Dw   =  ones(dimM*dimL, Ndat)*ew;
Dad  =  ones(dimP*dimL+dimP, Ndat)*ea;
DAll =  [Dd+Dad;Dw];
De   =  max(svd(DAll*DAll')) * eye(dimP*dimL+dimP+dimM*dimL);
dDe  =  ( (dimL+1) * ( max(svd(ones(dimP,1)*ed)) + max(svd(ones(dimP,1)*ea)) )^2 + ...
        dimL * max(svd(ones(dimM,1)*ew))^2 ) * Ndat;

tm1 =  dimP;
tm2  = dimP*dimL+dimP+dimM*dimL;
De11 = De(1:tm1,     1:tm1); De12 = De(1:tm1,     tm1+1:tm2); 
De21 = De(tm1+1:tm2, 1:tm1); De22 = De(tm1+1:tm2, tm1+1:tm2); 

Adz  = X0*X0'-De22;
Bdz  = -Lell'*X1*X0' + De12;
Cdz  = (Lell'*X1)*(Lell'*X1)' - De11;

eigDe22 = eig(De22);
eigAdz  = eig(Adz);



%% find the controller ====================================================
[vP, vK] = findControllerPetersenMIMO(dimL, dimP, dimM, Adz, Bdz, Cdz, Fell, Bell, Lell);
disp(vK);


% vK = [ -0.2896    0.2994   -0.1534   -0.3434    0.8220    0.7921   -0.5571   -0.2615;
%        -0.0129    0.0132   -0.1140    0.1189    0.0179    0.0172    0.1625   -0.0090 ];

eig_vK = eig(Aell+Bell*vK)

close_FLB = [Aini,       zeros(3,1), Bini*vK;
             zeros(1,3), Aa,         Ba*vK;
             Lell*Cini,  Lell*Ca,    Fell+Bell*vK];

eig_close_FLB = eig(close_FLB)




%% test1 the controller ===================================================
N = 81;
t = 0:1:N-1;

Uactual = zeros(dimM, N+dimL); % actual input
Unoise  = (rand(dimM, N+dimL)-0.5)*2* ew; % input noise
Umeasur = zeros(dimM, N+dimL); % measured input

Yartifi = zeros(dimP, N+dimL); % output of artificial system
Yactual = zeros(dimP, N+dimL); % actual output
Ynoise  = (rand(dimP, N+dimL)-0.5)*2* ed; % output noise
Ymeasur = zeros(dimP, N+dimL); % measured output

xoini = randn(dimN,1)*2;
xaini = zeros(dimNa,1);
xo = xoini; % actual system state with randn initial value
xa = xaini; % artificial system state with zero initial  

% implementation of the controller
% closed-loop system
for k = dimL+1 : N+dimL
    Yactual(:,k) = Cini*xo;
    Ymeasur(:,k) = Yactual(:,k) + Ynoise(:,k);
    
    Yartifi(:,k) = Ca*xa;

    Umeasur(:,k) = vK*[ Ymeasur(:,k-2)+Yartifi(:,k-2); 
                        Ymeasur(:,k-1)+Yartifi(:,k-1);
                        Umeasur(:,k-2);
                        Umeasur(:,k-1)                ];
    Uactual(:,k) =  Umeasur(:,k) - Unoise(:,k);

    xo = Aini*xo+Bini*Uactual(:,k);
    xa = Aa*xa+Ba*Umeasur(:,k);
end


figure('color', [1 1 1]);
set(gcf, 'Position', [100, 100, 500, 200]);

plot(t, Yactual(1,dimL+1:N+dimL), 'b-*'); hold on;
plot(t, Yactual(2,dimL+1:N+dimL), 'r--o'); hold on;
set(gca,'fontsize',12,'fontname','Times');
% set(gca,'XLim',[1 101]);
xlabel('steps','Fontname', 'Times New Roman','FontSize',12);
ylabel('output', 'Fontname', 'Times New Roman','FontSize',12);
legend('$y_1(k)$',...
       '$y_2(k)$',...
       'Fontname', 'Times New Roman','Interpreter', 'latex','FontSize',12);
grid on;



statf = 31;
figure('color', [1 1 1]);
set(gcf, 'Position', [100, 100, 300, 100]);

plot(t(statf:N), Yactual(1,dimL+statf:N+dimL), 'b-*'); hold on;
plot(t(statf:N), Yactual(2,dimL+statf:N+dimL), 'r--o'); hold on;
set(gca,'fontsize',12,'fontname','Times');
set(gca,'XLim',[t(statf) t(N)]);
grid on;

Ytemp = Yactual(:,dimL+1:N+dimL);


%% test2 the controller ===================================================
N = 81;
t = 0:1:N-1;

Unoise = Unoise(:,dimL+1:N+dimL);
Ynoise = Ynoise(:,dimL+1:N+dimL);

Uactual = zeros(dimM, N); % actual input
Umeasur = zeros(dimM, N); % measured input

Yartifi = zeros(dimP, N); % output of artificial system
Yactual = zeros(dimP, N); % actual output
Ymeasur = zeros(dimP, N); % measured output



xo = xoini; % actual system state with randn initial value
xa = xaini; % artificial system state with zero initial  
Ch = zeros(dimP*dimL+dimM*dimL, 1);

% implementation of the controller
% closed-loop system
for k = 1:N
    Yactual(:,k) = Cini*xo;
    Ymeasur(:,k) = Yactual(:,k) + Ynoise(:,k);
    
    Yartifi(:,k) = Ca*xa;

    Umeasur(:,k) = vK*Ch;
    Uactual(:,k) =  Umeasur(:,k) - Unoise(:,k);

    xo = Aini*xo+Bini*Uactual(:,k);
    xa = Aa*xa+Ba*Umeasur(:,k);
    Ch = Fell*Ch + Lell*( Ymeasur(:,k) + Yartifi(:,k) ) + Bell*Umeasur(:,k);
end

figure('color', [1 1 1]);
set(gcf, 'Position', [100, 100, 500, 200]);

plot(t, Yactual(1,1:N), 'b-*'); hold on;
plot(t, Yactual(2,1:N), 'r--o'); hold on;
set(gca,'fontsize',12,'fontname','Times');
% set(gca,'XLim',[1 101]);
xlabel('steps','Fontname', 'Times New Roman','FontSize',12);
ylabel('output', 'Fontname', 'Times New Roman','FontSize',12);
legend('$y_1(k)$',...
       '$y_2(k)$',...
       'Fontname', 'Times New Roman','Interpreter', 'latex','FontSize',12);
grid on;



statf = 31;
figure('color', [1 1 1]);
set(gcf, 'Position', [100, 100, 300, 100]);

plot(t(statf:N), Yactual(1,statf:N), 'b-*'); hold on;
plot(t(statf:N), Yactual(2,statf:N), 'r--o'); hold on;
set(gca,'fontsize',12,'fontname','Times');
set(gca,'XLim',[t(statf) t(N)]);
grid on;



































