function [go_r,t_r,s_r,go_d,t_d]=fslqr_reg(sys1,syscr,sysc1,sysc2)
%
% function [go_r,t_r,s_r,go_d,t_d] =  fslqr_reg(G,Cr,C1,C2)
% Frequency Shape LQR with reference input
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
%
% R       E                   U*    U        X         Y
% ---->0--->|Cr|-->O--->O------------->|G|------>|C|------->
%     -^          -^   -^         |           |        |
%      |           |    |         |           |        |
%      |           |     ---|C2|--            |        |
%      |           |                          |        |
%      |           -------------------|C1|-----        |
%      |-----------------------------------------------|
%
% Outputs
% go_r: open loop transfer function Go_r(s) with loop broken at E
% t_r: close-loop complementary sensitivity transfer function T_r(s)
% Y(s) = Tr(s) R(s)
% s_r: close-loop  sensitivity transfer function S(s)
% S_s(s) = I - T_r(s)
%
% go_d: open loop transfer function Go_d(s) with loop broken at input U
% -U -> U*
% t_d: close-loop complementary sensitivity transfer function T_d(s)
% T_d(s) = inv(I+Go_d(s))*Go_d(s)

sys1=minreal(sys1); syscr=minreal(syscr); sysc1=minreal(sysc1); 
sysc2=minreal(sysc2);
[a,b,c,d]=ssdata(sys1);
if norm(d)~=0
    error('plant must be strictly causal')
    return
end    

% open loop at the output
[n1,m1]=size(b);
[p1,n1]=size(c);
sys1s=ss(a,b,eye(size(a)),0*d);
sysc2c=feedback(eye(m1),sysc2);
sysgom=series(sysc2c,sys1s);
sysc1c=feedback(sysgom,sysc1);
sysgo=series(syscr,sysc1c);
sysgo=minreal(series(sysgo,c));
sysct=minreal(feedback(sysgo,eye(p1)));
syscs=minreal(feedback(eye(p1),sysgo));
go_r=sysgo;
t_r=sysct;
s_r=syscs;


% % open loop at control input

b1=minreal(series(c,syscr));
b2=minreal(parallel(b1,sysc1));
b3=minreal(series(sys1s,b2));
go_d=minreal(series(b3,sysc2c));
go_d=minreal(series(go_d,1*eye(m1)));
t_d = minreal(feedback(go_d,eye(m1)));

% b1=minreal(series(c,syscr));
% b2=minreal(parallel(b1,sysc1));
% b3=minreal(series(sys1s,b2));
% go_d=minreal(series(b3,sysc2c));
% t_d = minreal(feedback(go_d,eye(m1)));
% 


