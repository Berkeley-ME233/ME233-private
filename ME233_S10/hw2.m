clear
close all
clc
%1.a.
A = [1.2 0 0; 0 0 1; 0 -2 0];
B = [1 0; 0 0; 0 1];
C = [10 0 0];
R = [10 0; 0 10];
Q = C'*C;
Pinf = dare(A,B,Q,R)
N = 50;
x0 = [1 0 1]';
y0 = C*x0;
%i.
S = zeros(3);
P{N+1} = S;
J(N+1) = x0'*P{N+1}*x0/2;
for k = N+1:-1:2
    P{k-1} = Q+A'*P{k}*A-A'*P{k}*B*inv(R+B'*P{k}*B)*B'*P{k}*A;
    %K(k) = inv(R+B'*P(k)*B)*B'*P(k)*A;
    J(k-1) = x0'*P{k-1}*x0/2;
end
%P0 = Q+A'*P(1)*A-A'*P(1)*B*inv(R+B'*P(1)*B)*B'*P(1)*A;
%K0 = inv(R+B'*P0*B)*B'*P0*A;
%u0 = -K0*x0;
%x(1) = A*x0 + B*u0;
%u(1) = -K(2)*x(1);
figure
subplot(411)
plot(N:-1:0,J)
title('J^o[x_0,N-m,S,N]')
legend('S = 0','Location','Best')
P = P(1:2);
celldisp(P)
P{2}-P{1}
%ii.
S = diag([0 0 1]);
P{N+1} = S;
J(N+1) = x0'*P{N+1}*x0/2;
for k = N+1:-1:2
    P{k-1} = Q+A'*P{k}*A-A'*P{k}*B*inv(R+B'*P{k}*B)*B'*P{k}*A;
    J(k-1) = x0'*P{k-1}*x0/2;
end
subplot(412)
plot(N:-1:0,J)
legend('S = diag([0 0 1])','Location','Best')
P = P(1:2);
celldisp(P)
P{2}-P{1}
%iii.
S = diag([1 1 1]);
P{N+1} = S;
J(N+1) = x0'*P{N+1}*x0/2;
for k = N+1:-1:2
    P{k-1} = Q+A'*P{k}*A-A'*P{k}*B*inv(R+B'*P{k}*B)*B'*P{k}*A;
    J(k-1) = x0'*P{k-1}*x0/2;
end
subplot(413)
plot(N:-1:0,J)
legend('S = diag([1 1 1])','Location','Best')
P = P(1:2);
celldisp(P)
P{2}-P{1}
%iv.
S = diag([10 1 1]);
P{N+1} = S;
J(N+1) = x0'*P{N+1}*x0/2;
for k = N+1:-1:2
    P{k-1} = Q+A'*P{k}*A-A'*P{k}*B*inv(R+B'*P{k}*B)*B'*P{k}*A;
    J(k-1) = x0'*P{k-1}*x0/2;
end
subplot(414)
plot(N:-1:0,J)
xlabel('m')
legend('S = diag([10 1 1])','Location','Best')
P = P(1:2);
celldisp(P)
P{2}-P{1}
set(gcf,'Position',[360 60 560 670])

%1.b.
A = [1.2 1 0; 0 0 1; 0 -2 0];
Pinf = dare(A,B,Q,R)
%i.
S = zeros(3);
P{N+1} = S;
J(N+1) = x0'*P{N+1}*x0/2;
for k = N+1:-1:2
    P{k-1} = Q+A'*P{k}*A-A'*P{k}*B*inv(R+B'*P{k}*B)*B'*P{k}*A;
    J(k-1) = x0'*P{k-1}*x0/2;
end
figure
subplot(211)
plot(N:-1:0,J)
title('J^o[x_0,N-m,S,N]')
legend('S = 0','Location','Best')
P = P(1:2);
celldisp(P)
P{2}-P{1}
%ii.
S = diag([0 0 1]);
P{N+1} = S;
J(N+1) = x0'*P{N+1}*x0/2;
for k = N+1:-1:2
    P{k-1} = Q+A'*P{k}*A-A'*P{k}*B*inv(R+B'*P{k}*B)*B'*P{k}*A;
    J(k-1) = x0'*P{k-1}*x0/2;
end
subplot(212)
plot(N:-1:0,J)
xlabel('m')
legend('S = diag([0 0 1])','Location','Best')
P = P(1:2);
celldisp(P)
P{2}-P{1}

%1.c.
A = [0.8 1 1; 0 -1 1; 0 0 -1];
B = [0; 0; 1];
C = [0 1 0];
R = 0.1;
Q = C'*C;
Pinf = dare(A,B,Q,R)
%i.
S = zeros(3);
P{N+1} = S;
J(N+1) = x0'*P{N+1}*x0/2;
for k = N+1:-1:2
    P{k-1} = Q+A'*P{k}*A-A'*P{k}*B*inv(R+B'*P{k}*B)*B'*P{k}*A;
    J(k-1) = x0'*P{k-1}*x0/2;
end
figure
subplot(211)
plot(N:-1:0,J)
title('J^o[x_0,N-m,S,N]')
legend('S = 0','Location','Best')
P = P(1:2);
celldisp(P)
P{2}-P{1}
%ii.
S = diag([0 0 1]);
P{N+1} = S;
J(N+1) = x0'*P{N+1}*x0/2;
for k = N+1:-1:2
    P{k-1} = Q+A'*P{k}*A-A'*P{k}*B*inv(R+B'*P{k}*B)*B'*P{k}*A;
    J(k-1) = x0'*P{k-1}*x0/2;
end
subplot(212)
plot(N:-1:0,J)
xlabel('m')
legend('S = diag([0 0 1])','Location','Best')
P = P(1:2);
celldisp(P)
P{2}-P{1}

%3.a.
clear
%clc
close all
G = zpk([0 -2],[1 -.5 2],1,-1)
figure
rlocus(G'*G)
axis([-2.5 2.5 -2 2])
A = [2.5 -0.5 -1; 1 0 0; 0 1 0];
B = [1; 0; 0];
C = [1 2 0];
[P0 L0 K0] = dare(A,B,C'*C,32)
G0 = tf(ss(A,B,K0,0,-1))
figure
margin(G0)
figure
nyquist(G0)
axis([-2 0 -2 2])
figure
rlocus(G0)
axis([-2.5 2.5 -2 2])
title('Root locus of G_0(z)')
r = sqrt(32 / (32 + B'*P0*B))
2*asin(r/2)
2*asin(r/2)*180/pi