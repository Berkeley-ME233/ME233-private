clc
clear
close all

%%
%1.a)
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
%1.b)
G_a = tf(1,[1 0]);
G_e = G_p*G_a;
G_EA = G_PA*G_a;
%1.c)i.
figure
bode(G_e)
title('G_w(j \omega)')
set(gcf,'Position',[380 105 560 360])
%1.c)ii.
A_e = getfield(ss(G_e),'a');
B_e = getfield(ss(G_e),'b');
C_e = getfield(ss(G_e),'c');
figure
linestyles = {'b-','g--','r:'};
for mu=[1e-4 1e-3 1e-2]
    V = mu^2;
    M = care(A_e',C_e',B_e*B_e',V);
    L_e = M*C_e'*inv(V);
    G_okf = ss(A_e,L_e,C_e,0);
    mag = 20*log10(bode(G_okf,10^3));
    %allmargin(G_okf)
    bode(G_okf,linestyles{1})
    hold all
    linestyles(1) = [];
end
title('G_{okf}(j \omega)')
set(legend('\mu = 10^{-4}','\mu = 10^{-3}','\mu = 10^{-2}'),...
    'Position',[0.7 0.686 0.17 0.144])
%1.c)iii.
mu = 0.0064;
V = mu^2;
M = care(A_e',C_e',B_e*B_e',V);
L_e = M*C_e'*inv(V);
G_okf = ss(A_e,L_e,C_e,0);
figure
bodemag(1/DELTA,G_okf,'g--',G_okf*10^(5/20),'r:',{1e2,1e4})
legend('1/\Delta(j \omega)','G_{okf}(j \omega) with \mu = 0.0064',...
    'G_{okf}(j \omega) + 5 dB')
set(gcf,'Position',[380 105 560 260])
mag = 20*log10(bode(G_okf,10^3));
%1.c)v.
figure
step(feedback(G_okf,1))
set(gcf,'Position',[380 105 560 300])
set(title('Closed-Loop Step Response, step(feedback(G_okf,1))'),...
    'interpreter','none')
%%
%close all
%1.d)
rho = 1e-20;
[P eig_lqr K_e] = care(A_e,B_e,C_e'*C_e,rho);
K_e(1);
K_e(2);
K_e(3);
C_LQG = ss(A_e-B_e*K_e-L_e*C_e, L_e, K_e, 0);
G_o = G_e*C_LQG;
figure
margin(G_okf)
hold all
bode(G_o,'g--',{1e-1,1e3})
set(legend('G_{okf}(j \omega)','G_o(j \omega) for \rho = 10^{-20}'),...
    'Position',[0.605 0.7 0.273 0.133])
set(gcf,'Position',[380 105 560 360])
figure
step(feedback(G_o,1))
set(gcf,'Position',[380 105 560 300])
set(title('Closed-Loop Step Response, step(feedback(G_o,1))'),...
    'interpreter','none')
%1.e)
T_p = G_o/(1+G_o);
figure
bodemag(1/DELTA,T_p,'g--',{1e1,1e5})
legend('1/\Delta(j \omega)','T_p(j \omega)')
set(gcf,'Position',[380 105 560 360])
%%
close all
%1.f)i.
figure
margin(G_o)
figure
margin(G_PA*G_a*C_LQG)
hold all
bode(G_o,'g--',{1e0,1e5})
set(gcf,'Position',[380 105 560 360])
set(legend('G_{PA}(s)G_a(s)C_{LQG}(s)','G_p(s)G_a(s)C_{LQG}(s)'),...
    'Position',[0.226,0.194,0.273,0.128])
%1.f)ii.
figure
step(G_PA*G_a*C_LQG/(1+G_PA*G_a*C_LQG),G_o/(1+G_o),'g--')
set(gcf,'Position',[380 105 560 300])
axis([0 2 0 1.1])
set(legend('G_{PA}(s)G_a(s)C_{LQG}(s) / (1 + G_{PA}(s)G_a(s)C_{LQG}(s))',...
    'G_p(s)G_a(s)C_{LQG}(s) / (1 + G_p(s)G_a(s)C_{LQG}(s))'),...
    'Position',[0.364 0.426 0.505 0.177])
%1.g)i.
mu = 0.001
V = mu^2;
M = care(A_e',C_e',B_e*B_e',V);
L_e = M*C_e'*inv(V);
L_e(1)
L_e(2)
L_e(3)
G_okf = ss(A_e,L_e,C_e,0);
rho = 1e-15
[P eig_lqr K_e] = care(A_e,B_e,C_e'*C_e,rho);
K_e(1)
K_e(2)
K_e(3)
C_LQG = ss(A_e-B_e*K_e-L_e*C_e, L_e, K_e, 0);
G_o = G_e*C_LQG;
figure
step(G_PA*G_a*C_LQG/(1+G_PA*G_a*C_LQG))
title(['Step Response of G_{PA}(s)G_a(s)C_{LQG}(s) / '...
    '(1 + G_{PA}(s)G_a(s)C_{LQG}(s)) with \mu = 0.001, \rho = 10^{-15}'])
set(gcf,'Position',[380 105 560 300])
axis([0 2 0 1.1])
%1.g)ii.
mu = 1e-4
V = mu^2;
M = care(A_e',C_e',B_e*B_e',V);
L_e = M*C_e'*inv(V);
L_e(1)
L_e(2)
L_e(3)
G_okf = ss(A_e,L_e,C_e,0);
rho = 1e-14
[P eig_lqr K_e] = care(A_e,B_e,C_e'*C_e,rho);
K_e(1)
K_e(2)
K_e(3)
C_LQG = ss(A_e-B_e*K_e-L_e*C_e, L_e, K_e, 0);
G_o = G_e*C_LQG;
figure
step(G_PA*G_a*C_LQG/(1+G_PA*G_a*C_LQG))
title(['Step Response of G_{PA}(s)G_a(s)C_{LQG}(s) / '...
    '(1 + G_{PA}(s)G_a(s)C_{LQG}(s)) with \mu = 10^{-4}, \rho = 10^{-15}'])
set(gcf,'Position',[380 105 560 300])
axis([0 2 0 1.1])