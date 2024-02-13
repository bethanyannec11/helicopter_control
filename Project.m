clc, clear all, close all, format compact

% Rotorcraft dynamics
A = [-0.1778,zeros(1,4),-9.7807,-9.7807,zeros(1,4);
     0,-0.3104,0,0,9.7807,0,0,9.7807,zeros(1,3);
     -0.3326,-0.5353,zeros(1,4),75.7640,343.86,zeros(1,3);
     0.1903,-0.294,zeros(1,4),172.62,-59.958,zeros(1,3);
     0,0,1,zeros(1,8);zeros(1,3),1,zeros(1,7);
     zeros(1,3),-1,0,0,-8.1222,4.6535,zeros(1,3);
     0,0,-1,zeros(1,3),-0.0921,-8.1222,zeros(1,3);
     zeros(1,6),17.168,7.1018,-0.6821,-0.1070,0;
     0,0,-0.2834,zeros(1,5),-0.1446,-5.5561,-36.674;
     zeros(1,9),2.7492,-11.1120];
 
B = [zeros(6,3);
     0.0632,3.339,0;            % delta_lat (roll motions)
     3.1739,0.2216,0;           % delta_long (pitch motions)
     zeros(1,3);
     0,0,-74.364;               % delta_ped (yaw motions)
     zeros(1,3)];

C = [zeros(1,4),1,zeros(1,6);   % phi
     zeros(1,5),1,zeros(1,5);   % theta
     zeros(1,9),0,1];           % r_fb

D = zeros(size(C,1),size(B,2));

G = ss(A,B,C,D);

%% A
p = sigmaoptions("cstprefs") ;
p.MagUnits = 'abs' ;

So = minreal(inv(eye(3)+G)) ;
To = minreal(G*inv(eye(3)+G)) ;
Ti = minreal(G*inv(eye(3)+G)) ;
Si = minreal(inv(eye(3)+G)) ;


% Performance Functions
A = .005 ;
M = 2 ;
wb = 1 ;
wp = tf([1/M wb],[1 wb*A]) ;
Wp = wp*eye(3) ;

% Uncertainty Function
Wi = tf(makeweight(.1,50,2))*eye(3) ;
% Wi = tf([1 .2],[1 50])*eye(3) ;

% Nominal Plant
systemnames = 'G Wp' ;
inputvar = '[w(3); u(3)]' ;
outputvar = '[Wp; -G-w]' ;
input_to_G = '[u]' ;
input_to_Wp = '[w+G]' ;
sysoutname = 'P' ;
sysic ;
P = minreal(ss(P)) ;


%% B: Nominal Performance
n_meas = 3 ;
n_ctrl = 3 ;
w = logspace(-3,3) ;

[K,CL,gamma,info] = hinfsyn(P,n_meas,n_ctrl) ;
So = minreal(inv(eye(3)+G*K)) ;
To = minreal(G*K*inv(eye(3)+G*K)) ;
Ti = minreal(K*G*inv(eye(3)+K*G)) ;
Si = minreal(inv(eye(3)+K*G)) ;

GK = G*K ;


% Open Loop Sigma Values
figure(1)
sgtitle('Open Loop Singular Values')
subplot(2,1,1)
sigma(G(1,1),'b',G(2,2),'r--',p) 
hold on
title('\sigma(G)')
legend('\theta','\phi')

subplot(2,1,2)
sigma(GK(1,1),'b',GK(2,2),'r--',p) 
title('\sigma(GK)')
legend('\theta','\phi')

% CL TFs
figure(2)
subplot(2,1,1)
sigma(So(1,1),'b',So(2,2),'r--',inv(Wp),p)
legend('\theta','\phi','\sigma(1/Wp)')
title('Disturbance Rejection (S_o)')

WpSo = Wp*So ;
subplot(2,1,2)
sigma(WpSo(1,1),'b',WpSo(2,2),'r--',p)
legend('\theta','\phi')
title('Nominal Performance (NP)')

figure(3)
sigma(To(1,1),'b',To(2,2),'--r',p)
legend('\theta','\phi')
title('Tracking (T_o)')

figure(4)
WiTi = -Wi*Ti ;
sigma(WiTi(1,1),'b',WiTi(2,2),'--r',p)
legend('\theta','\phi')
title('Robust Stability (RS)')

% CL Step Response and Disturbance Rejection
figure(5)
subplot(2,2,1)
step(To(1,1))
ylabel('\theta')

subplot(2,2,2)
step(To(2,2))
ylabel('\phi')

subplot(2,2,3)
step(So(1,1))
ylabel('d_0 \rightarrow \theta')

subplot(2,2,4)
step(So(2,2))
ylabel('d_0 \rightarrow \phi')

%% C: Robust Performance
Wi = tf([1 .2],[.5 1])*eye(3) ;

A = .005 ;
M = 2 ;
wb = 1 ;
wp = tf([1/M wb],[1 wb*A]) ;
Wp = wp*eye(3) ;

% Nominal Plant
systemnames = 'G Wp Wi' ;
inputvar = '[ydel(3); w(3); u(3)]' ;
outputvar = '[Wi; Wp; -G-w]' ;
input_to_G = '[u+ydel]' ;
input_to_Wp = '[w+G]' ;
input_to_Wi = '[u]' ;
sysoutname = 'P' ;
sysic ;
P = minreal(ss(P)) ;

n_meas = 3 ;
n_ctrl = 3 ;
w = logspace(-3,3) ;


[K,CL,gamma,info] = hinfsyn(P,n_meas,n_ctrl) ;
So = minreal(inv(eye(3)+G*K)) ;
To = minreal(G*K*inv(eye(3)+G*K)) ;
Ti = minreal(K*G*inv(eye(3)+K*G)) ;
Si = minreal(inv(eye(3)+K*G)) ;

N = [-Wi*Ti Wi*K*So; -Wp*G*Si Wp*So] ;

N_frd = frd(N,w) ;
blk = [1 1; 1 1; 1 1; 3 3] ;
[mu,info] = mussv(N_frd,blk) ;

WiTi_frd = frd(N(1:3,1:3),w) ;
blk_WiTi = [1 1; 1 1; 1 1] ;
[mu_WiTi] = mussv(WiTi_frd,blk_WiTi) ;

WpSo_frd = frd(N(4:6,4:6),w) ;
blk_WpSo = [3 3] ;
[mu_WpSo] = mussv(WpSo_frd,blk_WpSo) ;


% RS and NP
figure(7)
subplot(2,1,1)
sigma(mu_WiTi(1,1),'b',mu_WiTi(1,2),'--r',p)
title('Robust Stability (RS)')

subplot(2,1,2)
sigma(mu_WpSo(1,1),'b',mu_WpSo(1,2),'--r',p)
title('Nominal Performance (NP)')


% Plot of mu bounds
figure(8)
sigma(mu(1,1),'b',mu(1,2),'--r',p)
ylabel('\mu upper/lower bounds (abs)')
title('Robust Performance (RP) \mu Plot')
legend('upper','lower')


% CL Step Response and Disturbance Rejection
figure(9)
subplot(2,2,1)
step(To(1,1))
ylabel('\theta')

subplot(2,2,2)
step(To(2,2))
ylabel('\phi')

subplot(2,2,3)
step(So(1,1))
ylabel('d_0 \rightarrow \theta')

subplot(2,2,4)
step(So(2,2))
ylabel('d_0 \rightarrow \phi')


%% D: D-K Iteration

% Initialize
D = append(1,1,1,tf(eye(3)),tf(eye(3))) ;

% Initial K-Step
[K,CL,gamma,info] = hinfsyn(D*P*inv(D),n_meas,n_ctrl) ;

% Initial D-Step
Nf = frd(lft(P,K),w) ;
[mu_Nf,mu_Info] = mussv(Nf,blk) ;
mu_RP(1) = norm(mu_Nf(1,1),inf,1e-6) ;

figure(11)
sigmaplot(mu_Nf(1,1),p)
ylabel('\mu upper bounds (abs)')
title('Manual DK Iteration')

% D-K Iteration Loop
n = 10 ;

for i = 2:n
    % Fit resulting D-scales
    [dsysl,dsysr] = mussvunwrap(mu_Info) ;
    dsysl = dsysl/dsysl(3,3) ;
    di = fitfrd(genphase(dsysl(1,1)),3) ;
    Di = append(di,di,di,tf(eye(3)),tf(eye(3))) ;

    % K-Step
    [Ki,CL,gamma(i),info] = hinfsyn(Di*P*inv(Di),n_meas,n_ctrl) ;

    % D-Step
    Nf = frd(lft(P,Ki),w) ;
    [mu_Nf,mu_Info] = mussv(Nf,blk) ;
    mu_RP(i) = norm(mu_Nf(1,1),inf,1e-6) ;
    if i == 6
        WiTi_frd = frd(Nf(1:3,1:3),w) ;
        blk_WiTi = [1 1; 1 1; 1 1] ;
        [mu_WiTi] = mussv(WiTi_frd,blk_WiTi) ;

        WpSo_frd = frd(Nf(4:6,4:6),w) ;
        blk_WpSo = [3 3] ;
        [mu_WpSo] = mussv(WpSo_frd,blk_WpSo) ;
        Kfinal = Ki ;

        figure(18)
        subplot(2,1,1)
        sigma(mu_WiTi(1,1),'b',mu_WiTi(1,2),'--r',p)
        title('Robust Stability (RS)')

        subplot(2,1,2)
        sigma(mu_WpSo(1,1),'b',mu_WpSo(1,2),'--r',p)
        title('Nominal Performance (NP)')

        figure(12)
        sigma(mu_Nf(1,1),p)
        ylabel('\mu upper bounds (abs)')
        title('Optimal Manual DK Iteration')
    end

    % Add to Plot
    figure(11)
    hold on
    sigma(mu_Nf(1,1),p)
    ylabel('\mu upper bounds (abs)')
    title('Manual DK Iterations')
    
end
legend('1','2','3','4','5','6','7','8','9','10')

% Peak of mu vs. iteration
figure(13)
plot(mu_RP)
hold on
plot(gamma)
xlabel('Iteration')
ylabel('\mu and \gamma')
title('Robust Performance Level and H_{\infty} Cost')
legend('Peak of \mu','\gamma')

% CL TFs
So1 = minreal(inv(eye(3)+G*Kfinal)) ;
To1 = minreal(G*Kfinal*inv(eye(3)+G*Kfinal)) ;
GK = G*Kfinal ;


% Closed Loop Step Response and Disturbance Rejection
figure(14)
subplot(2,2,1)
step(To(1,1),To1(1,1))
ylabel('\theta')
legend('Original','DK')

subplot(2,2,2)
step(To(2,2),To1(2,2))
ylabel('\phi')
legend('Original','DK')

subplot(2,2,3)
step(So(1,1),So1(1,1))
ylabel('d_0 \rightarrow \theta')
legend('Original','DK')

subplot(2,2,4)
step(So(2,2),So1(2,2))
ylabel('d_0 \rightarrow \phi')
legend('Original','DK')