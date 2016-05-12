function [ go_r, t_r, s_r, go_d, t_d]=fslqr_reg_est(sys1,syscr,sysc1,sysc2,kest)
%
% function [go_r,t_r,s_r,go_d,t_d] =  fslqr_reg_est(G,Cr,C1,C2,kest)
% Frequency Shape LQR with reference input
%   
% Only that observer state xhat(t) is used instead of actual state x(t)
% ***MUST use kest from function kalman ******
%
%
% Assembles FSLQ with reference input controller
% described in the block diagram  in page 38 in Lecture 14 Slides of ME233
%(Frequency Shape LQR)
% 
%
% Inputs:
% G = C(sI-A)^{-1}B : plant (must be strictly causal)
% Cr, C1, C2: Compensators Cr(s), C1(s) and C2(s) in page 38 in
% Lecture 15 Notes
% kest: Kalman filter estimator - MUST use kest from function kalman
%
% R       E                        U*   U       X         Y
% ---->0--->|Cr|-->O--->O----------------->|G|---->|C|------->
%    - ^          -^  - ^         |  |                  |
%      |           |    |         |  |                  |
%      |           |     ---|C2|--   |--->| kest|<------| 
%      |           |                          |         |
%      |           |                   Xhat   |         |
%      |           --------------|C1|---------|         |
%      |------------------------------------------------| 
%
% Outputs
% go_r: open loop transfer function Go_r(s) with loop broken at E
% t_r: close-loop complementary sensitivity transfer function T_r(s)
% Y(s) = Tr(s) R(s)
% s_r: close-loop  sensitivity transfer function S_r(s)
% S_r(s) = I - T_r(s)
%
% go_d: open loop transfer function Go_d(s) with loop broken at input U
% from -U -> U*
% t_d: close-loop complementary sensitivity transfer function T_d(s)
% T_d(s) = inv(I+Go_d(s))*Go_d(s)



sys1=minreal(sys1); syscr=minreal(syscr); sysc1=minreal(sysc1); 
sysc2=minreal(sysc2); kest=minreal(kest);

[a,b,c,d]=ssdata(sys1);
[n1,m1]=size(b);
[p1,n1]=size(c);
sysaug1=append(eye(m1),sys1);               % form augmented system to feed to kest
sysaug1uy=series(ones(2*m1,m1),sysaug1);    % input u1, outputs [u1;y1]
[aest,best,cest,dest]=ssdata(kest);
[p1est,n1est]=size(cest);                   % get #outputs and #states of kest
[n1est,m1est]=size(best);                   % get #inputs and #states of kest
if p1est-n1est ~= p1
        error('incompatible kest and plant')
    return
end
ckesty=zeros(n1est,p1+n1est);
ckesty(:,p1+1:p1+n1est)=eye(n1est);
kestny=series(kest,ckesty);                   % get only xh from kest
sysc1kestny=series(kestny,sysc1);             % C1(s) Kest
sysc2c=feedback(eye(m1),sysc2);               %(I+sysc2)^{-1}
sysgom=series(sysc2c,sysaug1uy);
sysc1cuy=feedback(sysgom,sysc1kestny);
cyonly=zeros(p1:p1+m1);
cyonly(:,m1+1:m1+p1)=eye(p1);
sysc1c=series(sysc1cuy,cyonly);
go_r=minreal(series(syscr,sysc1c));
t_r=minreal(feedback(go_r,ss(eye(p1))));
s_r=minreal(feedback(ss(eye(p1)),go_r));

yvectkest= m1+1:m1est ;                          % vector to indicate number of y inputs
yvectsyscr= 1:p1 ;
uvect= 1:m1 ;
b1 = parallel(sysc1kestny,syscr,yvectkest,yvectsyscr,uvect,uvect);
b2=minreal(series(b1,sysc2c));
b3 = minreal(feedback(b2,eye(m1),uvect,uvect));
yinputonly=eye(m1est);
yinputonly=yinputonly(:,m1+1:m1est);
b4=minreal(series(yinputonly,b3));
go_d =minreal(series(sys1,b4));
go_d=minreal(series(go_d, 1*eye(m1)));
t_d = minreal(feedback(go_d,ss(eye(m1))));

% open loop at control input

% % Get Kestny with only u as the input
% temp  = eye(m1est);
% tempu =temp(:,1:m1);
% sysc1kestnyUinputonly=minreal(series(tempu,sysc1kestny));   % C1(s)xh output, u input 
% 
% 
% tempy(1:m1,:)= temp(:, m1+1:m1est);
% sysc1kestnyYinputonly=minreal(series(tempu,sysc1kestny));   % C1(s)xh output, y input 
% 
% b1=minreal(parallel(sysc1kestnyUinputonly,sysc2));          % C1(s) kestUinput + C2(s)
% b2 = ss(eye(m1));                                               % eye(number of inputs)
% b3=minreal(feedback(b2,b1));                                % inv(I + C1(s) kestUinput + C2(s))
% 
% b4=minreal(parallel(sysc1kestnyYinputonly,syscr));          % C1(s) kestnoyY + Cr(s)
% LQG_FS = minreal(series(b4,b3));                            % Freq Shaped LQG
% 
% go_d = minreal(series(sys1,LQG_FS));                        % go_d(s) Open Loop TF from D
% t_d = minreal(feedback(go_d,ss(eye(m1))));                       % T_d(s)    comp


