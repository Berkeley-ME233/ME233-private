clear
clc
close all
Ts = 0.5;
qi = 1/tf('z',Ts);
qi.Variable = 'z^-1';
G = tf(1,[1 1 0]);
Gd = c2d(G,Ts);
Gd.Variable = 'z^-1';
A = tf(Gd.den{1},1,Ts,'Variable','z^-1');
B = tf(Gd.num{1},1,Ts,'Variable','z^-1')/qi;
b0 = B.num{1}(1);
pd1 = exp(Ts*[-1+i;-1-i]*sqrt(2)/2);
Am = 1-sum(pd1)*qi+prod(pd1)*qi^2;
bmo = 1-sum(pd1)+prod(pd1);
pd2 = exp(Ts*[-1+sqrt(3)*i;-1-sqrt(3)*i]);
Acp = 1-qi*sum(pd2)+qi^2*prod(pd2);
s0 = (Acp.num{1}(2)-Gd.den{1}(2))/b0;
s1 = (Acp.num{1}(3)-Gd.den{1}(3))/b0;
S = s0 + qi*s1;
Bs = B/b0;
Ad = 1;
R = Ad*Bs;
T = (Acp/qi)/b0;
ud = 1-2*qi^25+2*qi^50-2*qi^75;
d = 0.5*qi^40;
yd = qi*bmo*ud/Am;
r = yd*T;
y = qi*(b0*r+B*Ad*d)/Acp;
u = (r - S*y)/R;
%u2 = (A*r/Bs-qi*b0*S*d)/Acp;
figure
subplot(211)
step(ud,'--',yd,'g',y,'r')
legend('u_d(k)','y_d(k)','y(k)')
ylabel('')
title('')
xlim([0 60])
subplot(212)
step(u)
ylabel('u(k)')
title('')
%set(gcf,'Position',[380 105 560 360])
%%
Ad = 1-qi;
diop = Ad*A;
s0 = (Acp.num{1}(2)-diop.num{1}(2))/b0;
s1 = (Acp.num{1}(3)-diop.num{1}(3))/b0;
s2 = -diop.num{1}(4)/b0;
S = s0 + s1*qi + s2*qi^2;
R = Ad*Bs;
y = qi*(b0*r+B*Ad*d)/Acp;
u = (r - S*y)/R;
%u2 = (A*r/Bs-qi*b0*S*d)/Acp;
figure
subplot(211)
step(ud,'--',yd,'g',y,'r')
legend('u_d(k)','y_d(k)','y(k)')
ylabel('')
title('')
xlim([0 60])
subplot(212)
step(u)
ylabel('u(k)')
title('')
set(gcf,'Position',[380 105 560 360])
%%
close all
Bs = 1;
Bu = B;
[Rp, S] = bezout(diop.num{1}, Bu.num{1}, Acp.num{1});
Rp = tf(Rp,1,Ts,'Variable','z^-1');
S = tf(S,1,Ts,'Variable','z^-1');
R = Rp*Ad*Bs;
Ac = Acp*Bs;
T = Acp*qi'*Bu'/sum(Bu.num{1})^2;
r = yd*T;
y = qi*B*(r+R*d)/Ac;
%have to shift to calculate u due to acausality of feedforward
udelay = qi*(r - S*y)/R;
tend = 60;
[ud t1] = step(ud,tend);
t1 = [-5; t1(floor(1:.5:length(t1)+0.5))];
ud = [0; 0; ud(floor(1:.5:length(ud)))];
[yd t2] = step(yd,tend);
t2 = [-5; t2(floor(1:.5:length(t2)+0.5))];
yd = [0; 0; yd(floor(1:.5:length(yd)))];
[y t3] = step(y,tend);
t3 = [-5; t3(floor(1:.5:length(t3)+0.5))];
y = [0; 0; y(floor(1:.5:length(y)))];
[udelay t4] = step(udelay,tend+Ts); %go one step further for plot's sake
t4 = [-5; t4(floor(1:.5:length(t4)+0.5))-Ts]; %shifted one step
udelay = [0; 0; udelay(floor(1:.5:length(udelay)))];
figure
subplot(211)
plot(t1,ud,'--',t2,yd,t3,y)
legend('u_d(k)','y_d(k)','y(k)')
xlim([-5 60])
subplot(212)
plot(t4,udelay)
ylabel('u(k)')
xlabel('Time (sec)')
xlim([-5 60])
%set(gcf,'Position',[380 105 560 360])
