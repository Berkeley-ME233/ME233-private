clear
clc
close all
sys1 = tf([1 -0.3],[1 -0.5],-1);
N = 5000;
load('hw4.mat','w')
%w = randn(N,1);
[y k] = lsim(sys1,w,0:N-1);
[cov_ww lags_ww] = xcov(w,w,10);
[cov_wy lags_wy] = xcov(w,y,10);
[cov_yw lags_yw] = xcov(y,w,10);
[cov_yy lags_yy] = xcov(y,y,10);
figure
plot(lags_yy,(lags_yy>=0).*(0.34/0.75*0.5.^lags_yy + ...
    0.055/0.375*(lags_yy==0))+(lags_yy<=0)*0.34/0.75.*0.5.^(-lags_yy),'.-', ...
    lags_yy,cov_yy/N,'--')
set(title('$$\Lambda_{YY}(l)$$ vs $$l$$'),'interpreter','latex')
grid on
legend('predicted','simulated')
set(gcf,'Position',[476 322 560 290])

clear
A = [-0.08 -1; 0.7 0.1];
B = [0.34; 0.3];
C = [0 3];
C*inv(eye(2)-A)*B*10
for k=0:40
    my(k+1)=C*(eye(2)-A^k)*inv(eye(2)-A)*B*10;
    if k==0
        mx{k+1}=[0;0];
        lx{k+1}=0.1*eye(2);
    else
        mx{k+1}=A*mx{k}+B*10;
        lx{k+1}=A*lx{k}*A'+B*B';
    end
    my2(k+1)=C*mx{k+1};
    ly0(k+1)=C*lx{k+1}*C'+0.5;
    ly5(k+1)=C*(A^5)*lx{k+1}*C';
end
figure
plot(0:k,my,'.-')
grid on
set(title('$$m_y(k)$$ vs $$k$$'),'interpreter','latex')
%set(gcf,'Position',[476 459 560 153])
set(gcf,'Position',[476 322 560 290])
figure
subplot(211)
plot(0:k,ly0,'.-')
grid on
set(title('$$\Lambda_{YY}(k,0)$$ vs $$k$$'),'interpreter','latex')
subplot(212)
plot(0:k,ly5,'.-')
grid on
set(title('$$\Lambda_{YY}(k,5)$$ vs $$k$$'),'interpreter','latex')
lxbar=dlyap(A,B*B')
lybar=[];
lybar2=[];
for k=-10:10
    if k<0
        lybar=[lybar C*lxbar*(A^(-k))'*C'];
    elseif k==0
        lybar=[lybar C*lxbar*C'+0.5];
    else
        lybar=[lybar C*(A^k)*lxbar*C'];
    end
    lybar2=[lybar2 C*(A^abs(k))*lxbar*C'+0.5*(k==0)];
end
figure
plot(-10:10,lybar,'.-')
grid on
set(title('$$\bar{\Lambda}_{YY}(j)$$ vs $$j$$'),'interpreter','latex')
ylim([-3 4])
set(gcf,'Position',[360 381 560 317])

clc
close all
clear
sys2 = ss(zpk(-0.2,[-0.4 -0.8],1,-1))
s0 = sys2.C*dlyap(sys2.A,sys2.B*(sys2.B'))*(sys2.C)'
s1 = sys2.C*sys2.A*dlyap(sys2.A,sys2.B*(sys2.B'))*(sys2.C)'
s2 = sys2.C*sys2.A^2*dlyap(sys2.A,sys2.B*(sys2.B'))*(sys2.C)'
a = inv([s0 s1; s1 s0])*[s1;s2]
syt = s0 - a(1)*s1 - a(2)*s2
N = 1e7;
w = randn(N,1);
tic
y = lsim(sys2,w,0:N-1);
toc
yt = y;
tic
for i=1:N
    yt(i) = y(i)-a(1)*y(max(i-1,1))-a(2)*y(max(i-2,1));
end
toc
M = 1000;
s0sim = sum(y(M:N).^2)/(N-M)
sytsim = sum(yt(M:N).^2)/(N-M)