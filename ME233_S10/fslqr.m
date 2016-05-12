function [syscr,sysc1,sysc2,Ke,Le] =  fslqr(sys1,qr,qf,rf,rho)
%
% function [Cr,C1,C2,Le,Eig] =  fslqr(G,Qr,Qf,Rf,rho)
% Frequency Shape LQR with reference input
%   
% Obtains optimal controller for
%  J = int_-\infty^\inty { (Qr Y)^*(Qr Y) + (Qf X)^*(Qf X) + rho
% (Rf U)^*(Rf U)} dw
% 
%returns the compensators for the FSLQ with reference input control
% described Lecture 14 FSLQ Lecture Slides in ME233
% (Frequency Shape LQR)
%
% Inputs:
% G = C(sI-A)^{-1}B : plant (must be strictly causal)
% Qr : Output filter (must have at least one pole)
% Qf : State filter ( can be 0*eye(n)) where n is the number of states)
% Rf : Input filter
% rho>0 : input scalar weight
%
% Outputs
% Cr, C1, C2: Compensators Cr(s), C1(s) and C2(s)  
% Lecture 14 FSLQR Lecture Slides
% Le: extended state feedback gain
% Eig: vector of close loop poles
%
%
% R       E                   U*    U        X         Y
% ---->0--->|Cr|-->O--->O------------->|G|------>|C|------->
%    - ^          -^   -^         |           |        |
%      |           |    |         |           |        |
%      |           |     ---|C2|--            |        |
%      |           |                          |        |
%      |           -------------------|C1|-----        |
%      |-----------------------------------------------|
%

qr = minreal(ss(qr)); qf = minreal(ss(qf)); rf=minreal(ss(rf));
sys1=minreal(ss(sys1));
[a,b,c,d]=ssdata(sys1);
if norm(d)~=0
    error('plant must be strictly causal')
    return
end    
[ar,br,cr,dr]=ssdata(qr);
[nr,mr]=size(br);
[pr,nr]=size(cr);
qrr=series(ss(c),qr);
[arr,brr,crr,drr]=ssdata(qrr);
[a1,b1,c1,d1]=ssdata(qf);
[a2,b2,c2,d2]=ssdata(rf);
[n,m]=size(b);
[nrr,mrr]=size(brr);
[prr,nrr]=size(crr);
[n1,m1]=size(b1);
[p1,n1]=size(c1);
[n2,m2]=size(b2);
[p2,n2]=size(c2);
if m ~= m2
    error('inconsistent dimensions in rf')
    return
end
if m1 ~= n
    error('inconsistent dimensions in qf')
    return
end
if nrr == 0
    error('qr must be strictly causal')
    return
end
ae=zeros(n+nrr+n1+n2,n+nrr+n1+n2);
be=zeros(n+nrr+n1+n2,m);
ce=zeros(prr+p1+p2,n+nrr+n1+n2);
de=zeros(prr+p1+p2,m);
ae(1:n,1:n)=a;
be(1:n,1:m)=b;
ce(1:prr,1:n)=drr;
ce(prr+1:prr+p1,1:n)=d1;
de(prr+p1+1:prr+p1+p2,1:m) = sqrt(rho)*d2;
if nrr >= 1
%     ae(n+1:n+nrr,1:n)=br;
    ae(n+1:n+nrr,1:n+nrr)=[brr arr];
    ce(1:prr,n+1:n+nrr) = crr;
end
if n1 >= 1
    ae(n+nrr+1:n+nrr+n1,1:m1)=b1;
    ae(n+nrr+1:n+nrr+n1,n+nrr+1:n+nrr+n1)=a1;
    ce(prr+1:prr+p1,n+nrr+1:n+nrr+n1) = c1;
end
if n2 >= 1
    ae(n+nrr+n1+1:n+nrr+n1+n2,n+nrr+n1+1:n+nrr+n1+n2)=a2;
    be(n+nrr+n1+1:n+nrr+n1+n2,m) = b2;
    ce(prr+p1+1:prr+p1+p2,n+nrr+n1+1:n+nrr+n1+n2) = sqrt(rho)*c2;
end
qe=ce'*ce;
ne=ce'*de;
re=de'*de;
syse=ss(ae,be,ce,de);
[Ke,Pe,Le] = lqr(ae,be,qe,re,ne);
K=Ke(:,1:n);
if nrr >= 1
    Kr=Ke(:,n+1:n+nrr);
    [ar,br,cr,dr]=ssdata(qr);
    syscr=ss(ar,br,Kr,zeros(p2,mr));
else
    syscr=ss(zeros(p2,mr));
end
if n1 >= 1
    K1=Ke(:,n+nrr+1:n+nrr+n1);
    sysc1=parallel(ss(K),ss(a1,b1,K1,0*K));
else
    sysc1=ss(K);
end
if n2 >= 1
    K2=Ke(:,n+nrr+n1+1:n+nrr+n1+n2);
    sysc2=ss(a2,b2,K2,zeros(p2,m2));
else
    sysc2=ss(zeros(p2,m2));
end

syscr =minreal(syscr);
sysc1 = minreal(sysc1);
sysc2 = minreal(sysc2);

