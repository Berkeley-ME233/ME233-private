clear
clc
close all
A = [-0.08 -1; 0.7 0.1];
B = [0.34; 0.3];
C = [0 3];
W = 1;
X0 = [0.1 0; 0 0.1];
%loop for e
for V = [0.5 0.05 5]
    V
    %calculate for a
    M{1} = X0;
    for k=1:11
        M{k+1} = A*M{k}*A'+B*W*B'-A*M{k}*C'*inv(C*M{k}*C'+V)*C*M{k}*A';
        trM(k) = trace(M{k});
        Z{k} = M{k}-M{k}*C'*inv(C*M{k}*C'+V)*C*M{k};
        trZ(k) = trace(Z{k});
        L_pre(k) = C*M{k}*C'+V;
        L_post(k) = V*inv(C*M{k}*C'+V)*V;
    end
    trM_ss = trM(k)
    trZ_ss = trZ(k)
    L_pre_ss = L_pre(k)
    L_post_ss = L_post(k)
    %plot for b
    figure
    subplot(4,1,1)
    plot(0:k-1,L_pre,'.-')
    grid on
    set(title(['$$V=' num2str(V) '$$']),'interpreter','latex')
    set(ylabel('$$\Lambda_{{\tilde{y}}^o {\tilde{y}}^o}(k,0)$$'),'interpreter','latex')
    subplot(4,1,2)
    plot(0:k-1,L_post,'.-')
    grid on
    set(ylabel('$$\Lambda_{\tilde{y} \tilde{y}}(k,0)$$'),'interpreter','latex')
    subplot(4,1,3)
    plot(0:k-1,trM,'.-')
    grid on
    set(ylabel('trace$$\lbrace M(k) \rbrace$$'),'interpreter','latex')
    subplot(4,1,4)
    plot(0:k-1,trZ,'.-')
    grid on
    set(xlabel('$$k$$'),'interpreter','latex')
    set(ylabel('trace$$\lbrace Z(k) \rbrace$$'),'interpreter','latex')
    %dare for c
    Mbar = dare(A',C',B*W*B',V)
    Zbar = Mbar-Mbar*C'*inv(C*Mbar*C'+V)*C*Mbar
    Lbar_pre = C*Mbar*C'+V
    Lbar_post = V*inv(C*Mbar*C'+V)*V
    Fbar = Mbar*C'*inv(C*Mbar*C'+V)
    Lbar = A*Fbar
    cl_eig = eig(A-Lbar*C)
    %simulation loop for d
    N = 10000;
    Nss = 1000;
    u = 10;
    x{1} = sqrt(0.1)*randn(2,1);
    xhato{1} = x{1};
    w = randn(N,1);
    v = sqrt(V)*randn(N,1);
    for k=1:N
        x{k+1} = A*x{k}+B*u+B*w(k);
        y(k) = C*x{k}+v(k);
        ytildeo(k) = y(k)-C*xhato{k};
        xhat{k} = xhato{k}+M{k}*C'*inv(C*M{k}*C'+V)*ytildeo(k);
        ytilde(k) = y(k)-C*xhat{k};
        M{k+1} = A*M{k}*A'+B*W*B'-A*M{k}*C'*inv(C*M{k}*C'+V)*C*M{k}*A';
        xhato{k+1} = A*xhato{k}+B*u+A*M{k}*C'*inv(C*M{k}*C'+V)*ytildeo(k);
        xtildeo{k} = x{k}-xhato{k};
        xtilde{k} = x{k}-xhat{k};
    end
    %covariances for e
    Mbarsim = cov([xtildeo{Nss:N}]')
    Zbarsim = cov([xtilde{Nss:N}]')
    Lbar_pre_sim = cov(ytildeo(Nss:N))
    Lbar_post_sim = cov(ytilde(Nss:N))
end