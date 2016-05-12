clear
close all
sys1 = tf([1 -0.3],[1 -0.5],-1);
N = 5000;
load('hw4.mat','w')
%w = randn(N,1);
[y k] = lsim(sys1,w,0:N-1);
[cov_ww lags_ww] = xcov(w,w,10,'coeff');
[cov_wy lags_wy] = xcov(w,y,10,'coeff');
[cov_yw lags_yw] = xcov(y,w,10,'coeff');
[cov_yy lags_yy] = xcov(y,y,10,'coeff');
figure
subplot(221)
plot(lags_ww,cov_ww,'.-')
title('\Lambda_{WW}(j)')
xlabel('j')
grid on
subplot(222)
plot(lags_wy,cov_wy,'.-')
title('\Lambda_{WY}(j)')
xlabel('j')
grid on
subplot(223)
plot(lags_yw,cov_yw,'.-')
title('\Lambda_{YW}(j)')
xlabel('j')
grid on
subplot(224)
plot(lags_yy,cov_yy,'.-')
title('\Lambda_{YY}(j)')
xlabel('j')
grid on
figure
subplot(211)
plot(lags_yy,(lags_yy>=0).*(0.4*0.5.^lags_yy+0.6*(lags_yy==0)),'.-')
title('\Lambda_{YW}(j) predicted')
xlabel('j')
grid on
subplot(212)
plot(lags_yy,(lags_yy<=0).*(0.4*0.5.^(-lags_yy)+0.6*(lags_yy==0)),'.-')
title('\Lambda_{WY}(j) predicted')
xlabel('j')
grid on
figure
subplot(211)
plot(lags_yy,(lags_yy>=0).*(0.4*0.5.^lags_yy+0.6*(lags_yy==0)),...
    lags_yy,cov_yw,'--')
title('\Lambda_{YW}(j)')
xlabel('j')
grid on
legend('predicted','simulated')
subplot(212)
plot(lags_yy,(lags_yy<=0).*(0.4*0.5.^(-lags_yy)+0.6*(lags_yy==0)),...
    lags_yy,cov_wy,'--')
title('\Lambda_{WY}(j)')
xlabel('j')
grid on
legend('predicted','simulated')