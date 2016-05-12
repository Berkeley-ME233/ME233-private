function [go_r,t_r,s_r,go_d,t_d]=fslqr_reg_robust_test(sys_act,sys_nom,syscr,sysc1,sysc2)
%
% function [go_r,t_r,s_r,go_d,t_d] =  fslqr_reg_robust_test(GA,G,Cr,C1,C2)
% Frequency Shape LQR with reference input
%   
% Assembles FSLQ with reference input controller
% described in the block diagram  in page 38 in Lecture 14 Slides of ME233
%(Frequency Shape LQR)
%
% and checks its robustness to multiplicative input uncertainties
% 
%
% Inputs:
% GA: Actual plant: plant (must be strictly causal)
% G = C(sI-A)^{-1}B : nominal plant (must be strictly causal)
% Cr, C1, C2: Compensators Cr(s), C1(s) and C2(s) in page 38 in
% Lecture 15 Notes
%
%  R       E                   U*     U                 X          Y
%  ---->0--->|Cr|-->O--->O------------->|I + D|-->|G|------>|C|------->
%      -^          -^   -^         |                    |        |
%       |           |    |         |                    |        |
%       |           |     ---|C2|--                     |        |
%       |           |                                   |        |
%       |           -------------------|C1|-------------|        |
%       |--------------------------------------------------------| 
%
% where D = inv(G)*[GA-G]
% 
% Outputs
% go_r: open loop transfer function Go_r(s) with loop broken at E
% t_r: close-loop complementary sensitivity transfer function T_r(s)
% Y(s) = Tr(s) R(s)
% s_r: close-loop  sensitivity transfer function S(s)
% S_s(s) = I - T_r(s)
%
% go_d: open loop transfer function Go_d(s) with loop broken at input U
% from -U -> U*
% t_d: close-loop complementary sensitivity transfer function T_d(s)
% T_d(s) = inv(I+Go_d(s))*Go_d(s)





I_plus_delta =  minreal(inv(sys_nom)*sys_act);
 

sys1=minreal(sys_nom); syscr=minreal(syscr); sysc1=minreal(sysc1); 
sysc2=minreal(sysc2);
[a,b,c,d]=ssdata(sys1);
if norm(d)~=0
    error('plant must be strictly causal')
    return
end    

% open loop at the output
[n1,m1]=size(b);
[p1,n1]=size(c);
if m1 ~= p1
        error('plant must be square')
    return
end

sys1s=ss(a,b,eye(size(a)),0*d);
sys1s=minreal(series(sys1s,I_plus_delta));
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


