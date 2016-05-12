function hw10p1
clc
close all

%%
%1.a)i.
zeta_b = 0.707;
omega_b = 10;
G_p = tf(omega_b^2,[1 2*zeta_b*omega_b omega_b^2]);
zeta_r = 0.015;
omega_r = 1000;
zeta_t = 0.015;
omega_t = 1200;
omega_n = 0.9*omega_t;
G_PA = G_p*omega_r*omega_t^2/omega_n^2*...
    tf([zeta_r omega_r],[1 2*zeta_r*omega_r omega_r^2])*...
    tf([1 2*zeta_t*omega_n omega_n^2],[1 2*zeta_t*omega_t omega_t^2]);
DELTA = G_PA/G_p - 1;
Q_r = tf(1,[1 0]);
[syscr,sysc1,sysc2,Ke,Le] = fslqr(G_p,Q_r,zeros(2),1,3.2e-9);
[go_r,t_r,s_r,go_d,t_d] = fslqr_reg(G_p,syscr,sysc1,sysc2);
figure
margin(go_r)
%1.a)ii.
figure
step(t_r)
title('Closed-loop nominal step response for \rho = 3.2e-9')
set(gcf,'Position',[380 105 560 200])
%1.a)iii.
figure
bodemag(go_d,t_d,'g--',1/DELTA,'r:')
title('Time (sec)')
axis([1 1e4 -50 100])
set(gcf,'Position',[380 105 560 360])
legend('G_{o(-A)\rightarrowB}(j \omega)',...
    'T_{(-A)\rightarrowB}(j \omega)','1/\Delta(j \omega)')
%1.a)iv.
[go_ra,t_ra,s_ra,go_da,t_da] = fslqr_reg_robust_test(G_PA,G_p,syscr,sysc1,sysc2);
clpoles = pole(t_ra);
%clpoles(1)
%clpoles(3)
%1.b)i.
close all
%coefs = fmincon(@(coefs) slowesteigenvalue(coefs,G_PA,G_p,Q_r),[1000 500 0.7 0.7],...
%    [],[],[],[],[0 0 0.25 0.25],[1000 1000 10 10],@(coefs) nonlcon(coefs,60,G_p,Q_r));
%coefs = [1000 250 1.6 0.25];
%p = coefs(1);
%z = coefs(2)
%xip = coefs(3)
%xiz = coefs(4)
%R_f = p^2/z^2*tf([1 2*xiz*z z^2],[1 2*xip*p p^2]);
R_f = 16*tf([1 125 250^2],[1 1414 1000^2]);
%R_f = 100/9*tf([1 60 300^2],[1 1600 1000^2]);
%R_f = 4*tf([1 700 500^2],[1 600 1000^2]);
%R_f = 100/49*tf([1 2100 700^2],[1 500 1000^2]);
%R_f = 4*tf([1 700 500^2],[1 500 1000^2]);
[syscr,sysc1,sysc2,Ke,Le] = fslqr(G_p,Q_r,zeros(2),R_f,3.2e-9);
[go_r,t_r,s_r,go_d,t_d] = fslqr_reg(G_p,syscr,sysc1,sysc2);
figure
margin(go_r)
figure
bodemag(go_d,t_d,'g--',1/DELTA,'r:')
title('Small gain condition with frequency weight R_f')
axis([1 1e4 -50 100])
set(gcf,'Position',[380 105 560 360])
legend('G_{o(-A)\rightarrowB}(j \omega)',...
    'T_{(-A)\rightarrowB}(j \omega)','1/\Delta(j \omega)')
%1.b)ii.
[go_ra,t_ra,s_ra,go_da,t_da] = fslqr_reg_robust_test(G_PA,G_p,syscr,sysc1,sysc2);
clpoles = pole(t_ra)
clpoles(2)
clpoles(4)
%1.c)i.
mu = 1e-8;
[A B C D] = ssdata(G_p);
kest = kalman(ss(A,[B B],C,[D D]),1,mu^2);
[go_rkf, t_rkf, s_rkf, go_dkf, t_dkf] = fslqr_reg_est(G_p,syscr,sysc1,sysc2,kest);
figure
bode(go_dkf,go_d,'g--')
title(['Loop Transfer Recovery for \mu = ' strrep(num2str(mu),'0','')])
set(gcf,'Position',[380 105 560 360])
set(legend('Kalman Filter G_{o(-A)\rightarrowB}(j \omega)',...
    'state feedback G_{o(-A)\rightarrowB}(j \omega)'),'Position',[0.18 0.51 0.34 0.13])
%1.c)ii.
figure
bodemag(go_dkf,t_dkf,'g--',1/DELTA,'r:')
title('Small gain condition with Kalman Filter')
axis([1 1e6 -150 100])
set(gcf,'Position',[380 105 560 360])
legend('G_{o(-A)\rightarrowB}(j \omega)',...
    'T_{(-A)\rightarrowB}(j \omega)','1/\Delta(j \omega)')
figure
margin(go_rkf)
hold all
bode(go_r,{1,1e8})
%1.c)iii.
[go_ra, t_ra, s_ra, go_da, t_da] = fslqr_reg_est(G_PA,syscr,sysc1,sysc2,kest);
figure
margin(go_ra)
xlim([10 1e6])
set(gcf,'Position',[380 105 560 360])
title(['Open-loop G_{o E\rightarrowY}(j \omega) with actual plant model' 10 ...
    'Gm = 9.02 dB (at 158 rad/sec) ,  Pm = 60.3 deg (at 55.3 rad/sec)'])
allmargin(go_ra)
figure
step(t_ra)
set(gcf,'Position',[380 105 560 300])
title('Closed-Loop Step Response of system with Kalman Filter and actual model')

%%
%{

function [c ceq] = nonlcon(coefs,wgc,G_p,Q_r)
p = coefs(1);
z = coefs(2);
xip = coefs(3);
xiz = coefs(4);
R_f = p^2/z^2*tf([1 2*xiz*z z^2],[1 2*xip*p p^2]);
[syscr,sysc1,sysc2,Ke,Le] = fslqr(G_p,Q_r,zeros(2),R_f,3.2e-9);
[go_r,t_r,s_r,go_d,t_d] = fslqr_reg(G_p,syscr,sysc1,sysc2);
[gm pm wp wg] = margin(go_r);
c = 0;
ceq = wg - wgc;

function lambda = slowesteigenvalue(coefs,G_PA,G_p,Q_r)
p = coefs(1);
z = coefs(2);
xip = coefs(3);
xiz = coefs(4);
R_f = p^2/z^2*tf([1 2*xiz*z z^2],[1 2*xip*p p^2]);
[syscr,sysc1,sysc2,Ke,Le] = fslqr(G_p,Q_r,zeros(2),R_f,3.2e-9);
[go_ra,t_ra,s_ra,go_da,t_da] = fslqr_reg_robust_test(G_PA,G_p,syscr,sysc1,sysc2);
lambda = max(real(pole(t_ra)));
disp(coefs)
disp(lambda)

%}