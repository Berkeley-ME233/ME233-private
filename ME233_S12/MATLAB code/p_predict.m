% ME 233, Spring 2006, sp_predict.m
% Parameter Estimation using Recursive least squares
% Series parallel predictor

% revised for Spring 2012 by Richard Conway

clear
% actual values of the parameters
disp('**** SERIES PARALLEL PREDICTOR ****');
disp('         b1*z+b2');
disp('   G(z)=-----------');
disp('        z^2+a1*z+a2');
disp('input the actual parameter values and hit return');

%   Maximum number of iterations
N = 3000;
t=0:N-1;
%   Plant parameters
a1=1.7;
a2=0.72;
b1=.1;
b2=0.05;

%Parallel estimator parameters
c1=1.3;
c2=0.42;

disp('');
disp('**** select the input function ****');
disp('if u=random, then type -1');
disp('if u=sin(w1*k)+sin(w2*k)+sin(w3*k), then type 0');
disp('if u=1, then type 1');
ct=input('  :  ');

if ct==0,
  disp('**** input values of w1, w2, w3 ****');
  w1=input('w1=');
  w2=input('w2=');
  w3=input('w3=');
  u = sin(w1*t)+sin(w2*t)+sin(w3*t);
elseif ct==-1
    u=rand(size(t));
elseif ct==1
    u=ones(size(t));
end

disp('');
disp('**** select the measurement noise function ****');
disp('if v= white noise, then type -1');
disp('if v=0, then type 0');
disp('if v= colored noise, then type 1');
ctv=input('  :  ');

if ctv==0,
    vvn=zeros(size(t));
    av=0;
    bv=1;
elseif ctv==1
    vvn=randn(size(t));
    av=.8;
    bv=.2;
elseif ctv==-1
    vvn=randn(size(t));
    av=0;
    bv=1;  
end

ct1=1;
while ct1==1,
    disp('**** input forgetting factor ****');
    disp('**** lamda1=1, lamda2=0: constant');
    disp('**** lamda1=1, lamda2=1: least squares');
    disp('**** lamda1<1, lamda2=1: weighted least squares');
    lam1=input('lamda1=');
    lam2=input('lamda2=');

    % initial gain matrix
    disp('**** input the diagonal elements of the gain matrix ****');
    disp('**** F=diag(f11, f22, f33, f44) ****');
    f11=input('f11=');
    f22=input('f22=');
    f33=input('f33=');
    f44=input('f44=');

    F=diag([f11, f22, f33, f44]);
    F2 = F;
    
    %preallocates arrays, if necessary
    if ~exist('thetahat','var')
        thetahat = zeros(N,4);
        yv = zeros(1,N);
        yhatv = zeros(1,N);
        err = zeros(1,N);
        f = zeros(1,N);
        t = 1:N;
    end

    % estimates of the parameters
    thhat=[0 0 0 0]';

    % initial data values
    % y1=y(k-1), y2=y(k-2), u1=u(k-1), u2=u(k-2)
    y1=0;
    y2=0;
    u1=0;
    u2=0;
    yh1=0;
    yh2=0;
    e1=0;
    e2=0;
    vn1=0;
    phihat=[0 0 0 0]';

    for k=1:N,
        % generate noise
        vn = av*vn1+bv*vvn(k);
        % make new measurements

        y=-a1*y1-a2*y2+b1*u1+b2*u2 + b2*vn;

        % find a-priori error
        yho=thhat'*phihat;
        eo=y-yho;
        vo=eo+c1*e1+c2*e2;
        
        % update parameter values
        Fp=F*phihat;
        pFp=phihat'*Fp;
        v = lam1/(lam1+pFp)*vo;
        thhat=thhat+Fp*v/lam1;

        % find a-posteriori error
        yh=thhat'*phihat;
        e=y-yh;

        % update the gain matrix
        F = (F - lam2*(Fp*Fp')/(lam1+lam2*pFp)) / lam1;
        
        % prepare for next iteration
        yh2=yh1;
        yh1=yh;
        e2=e1;
        e1=e;
        y2=y1;
        y1=y;
        u2=u1;
        u1=u(k);
        phihat=[-yh1 -yh2 u1 u2]';
        vn1=vn;

        % store values for plotting
        thetahat(k,1)=thhat(1); 
        thetahat(k,2)=thhat(2);
        thetahat(k,3)=thhat(3); 
        thetahat(k,4)=thhat(4);

        yv(k)=y;
        yhatv(k)=yho;
        err(k)=eo;
        f(k)=det(F);

    end

    % plotting error, gain matrix and estimates
    figure(1);
    subplot(211)
    plot(t,err); 
    title('prediction error, parallel predictor');
    subplot(212)
    plot(t, f);
    title('determinant of gain matrix');
    if ct==0,
        xlabel('input: sum of sinusoids');
    elseif ct==1,
        xlabel('input: constant');
    elseif ct==-1,
        xlabel('input: random');
    end % ending if

    %   disp('pause ...'); pause
    figure(2);
    subplot(221); 
    plot(t, thetahat(:, 1),'b', t , a1*ones(size(thetahat(:, 1))),'-k'); 
    ylabel('a1');
    subplot(222); 
    plot(t, thetahat(:, 2),'b', t , a2*ones(size(thetahat(:, 2))),'-k'); 
    ylabel('a2');
    subplot(223); 
    plot(t, thetahat(:, 3),'b', t , b1*ones(size(thetahat(:, 3))),'-k'); 
    ylabel('b1');
    subplot(224); 
    plot(t, thetahat(:, 4),'b', t , b2*ones(size(thetahat(:, 4))),'-k'); 
    ylabel('b2');

%     if ctv==0,
%         title('prediction error, parallel predictor, no noise ');
%     elseif ctv==1,
%         title('prediction error, parallel predictor, colored noise');
%     elseif ctv==-1,
%         title('prediction error, parallel predictor, white noise');
%     end
    ct1=input('do you wish to try different lambda values?, yes - type 1 else 0: ');
end