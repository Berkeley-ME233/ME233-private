close all; clear all; clc

load me233hw6p1_model

% alpha = [1 1 1];
alpha = [2 2 2];
% alpha = [1e-5 0.07 0.05];

%lqg cost
Q = C1'*C1 + alpha(1)*C2'*C2;
R = diag(alpha(2:3));

%LQR design
sys_lqr = ss(A,B,[],[],-1);
K = lqr(sys_lqr,Q,R);

%Kalman filter design
sys_kf = ss(A,[B Bw],[C1; C2],zeros(2,4),-1);
Kest = kalman(sys_kf, eye(2), eye(2), 'current');

%optimal LQG controller
Klqg = lqgreg(Kest,K);

%closed-loop system
Glft = ss(A, [Bw zeros(11,2) B], [C1; C2; zeros(2,11); C1; C2], ...
    [zeros(2,6); zeros(2,4) eye(2); zeros(2,2) eye(2) zeros(2,2)], -1);
Gcl = lft(Glft,Klqg);

%results
disp(['3 sigma PES: ' num2str(3*norm(Gcl(1,:)))])
disp(['3 sigma RPES: ' num2str(3*norm(Gcl(2,:)))])
disp(['3 sigma Uv: ' num2str(3*norm(Gcl(3,:)))])
disp(['3 sigma Um: ' num2str(3*norm(Gcl(4,:)))])


