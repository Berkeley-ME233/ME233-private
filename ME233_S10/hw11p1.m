clear
clc
close all
Ts = 1;
qi = tf([0 1],1,Ts,'Variable','z^-1');
A = (1-0.8*qi)*(1-0.7*qi);
d = 1;
B = 0.1;
Bs = 1;
N = 8;
Ad = 1-qi^N;
Ac = 1;
Rp = 1;
R = Rp*Ad*Bs;
S = (Ac - Ad*A*Rp)/(qi^d*B);
% disturbance input is the impulse response of the following:
dk = tf([0 1.5 3 0 -2 -2 0 0.5],1,Ts,'Variable','z^-1')/Ad;
%% 1.b
y = dk*feedback(qi^d*B/A,S/R);
u = -y*S/R;
tend = 25;
figure
subplot(311)
%impulse(dk,tend)
plot([floor(0.5+(0:0.5:tend)) tend],reshape(repmat(impulse(dk,tend)',[2 1]),1,[]))
ylabel('d(k)')
ylim([-3 4])
subplot(312)
%impulse(y,tend)
plot([floor(0.5+(0:0.5:tend)) tend],reshape(repmat(impulse(y,tend)',[2 1]),1,[]))
ylabel('y(k)')
ylim([-0.3 0.4])
subplot(313)
%impulse(u,tend)
plot([floor(0.5+(0:0.5:tend)) tend],reshape(repmat(impulse(u,tend)',[2 1]),1,[]))
xlabel('k')
ylabel('u(k)')
ylim([-4 4])
set(gcf,'Position',[380 105 560 400])
%% 1.c
linestyles = {'-','--',':'};
figure
for kr = [0.5 1 1.3]
    y = dk*feedback(qi^d*B/A,kr*qi^(N-d)*A/(Ad*B));
    u = -y*kr*qi^(N-d)*A/(Ad*B);
    tend = 50;
    subplot(3,1,1)
    %impulse(dk,tend)
    plot([floor(0.5+(0:0.5:tend)) tend], ...
        reshape(repmat(impulse(dk,tend)',[2 1]),1,[]),linestyles{1})
    ylabel('d(k)')
    ylim([-3 4])
    hold all
    subplot(3,1,2)
    %impulse(y,tend)
    plot([floor(0.5+(0:0.5:tend)) tend], ...
        reshape(repmat(impulse(y,tend)',[2 1]),1,[]),linestyles{1})
    ylabel('y(k)')
    ylim([-0.4 0.8])
    hold all
    subplot(3,1,3)
    %impulse(u,tend)
    plot([floor(0.5+(0:0.5:tend)) tend], ...
        reshape(repmat(impulse(u,tend)',[2 1]),1,[]),linestyles{1})
    xlabel('k')
    ylabel('u(k)')
    ylim([-5 4])
    hold all
    set(gcf,'Position',[380 105 560 500])
    linestyles(1) = [];
end
subplot(3,1,2)
legend('k_r = 0.5','k_r = 1','k_r = 1.3','Location','NE', ...
    'Orientation','Horizontal')
%% 1.d
Abar = (1-0.8*qi)*(1-0.7*qi);
G = 0.1*qi/Abar;
G_A = G*0.8*qi/(1-0.2*qi);
figure
rlocus(G_A*qi^(N-d)*Abar/(0.1*Ad))
axis([-1.2 1.2 -1.2 1.2])
axis equal
%% 1.e.i
Q = (qi^-1 + 2 + qi)/4;
figure
rlocus(G_A*qi^(N-d)*Abar/(0.1*(1-Q*qi^N)))
title('Root Locus with Q-filter')
axis([-1.2 1.2 -1.2 1.2])
axis equal
%% 1.e.ii
%tend = 250;
%obj = @(kr) sum(impulse(dk*feedback(G_A,kr*qi^(N-d)*Abar/(0.1*(1-Q*qi^N)))).^2);
obj = @(kr) max(abs(eig(feedback(kr*G_A*qi^(N-d)*Abar/(0.1*(1-Q*qi^N)),1))));
figure
plot(eps:1e-3:1,arrayfun(obj,eps:1e-3:1))
title('Max closed-loop eigenvalue magnitude vs k_r')
%ylabel('max(abs(\lambda))')
ylim([0.95 1])
set(gcf,'Position',[380 105 560 200])
kr = fminsearch(obj,1);
y = dk*feedback(G_A,kr*qi^(N-d)*Abar/(0.1*(1-Q*qi^N)));
u = -y*kr*qi^(N-d)*Abar/(0.1*(1-Q*qi^N));
tend = 100;
figure
subplot(311)
%impulse(dk,tend)
plot([floor(0.5+(0:0.5:tend)) tend],reshape(repmat(impulse(dk,tend)',[2 1]),1,[]))
title('Control with Q-filter, k_r = 0.444')
ylabel('d(k)')
ylim([-3 4])
subplot(312)
%impulse(y,tend)
plot([floor(0.5+(0:0.5:tend)) tend],reshape(repmat(impulse(y,tend)',[2 1]),1,[]))
ylabel('y(k)')
ylim([-0.25 0.75])
subplot(313)
%impulse(u,tend)
plot([floor(0.5+(0:0.5:tend)) tend],reshape(repmat(impulse(u,tend)',[2 1]),1,[]))
xlabel('k')
ylabel('u(k)')
ylim([-3 3])
set(gcf,'Position',[380 105 560 450])
