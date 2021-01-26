function [] = Compute_SensitivityActivation(u_v,p,xl,Col)
k = p(1);
th = p(2);
ut = p(7);
u_v = u_v*ut;
I = p(9);
% Return values of titration of A, I and its hill function coefficient.
A_v = A_F_u_I(u_v,p);
A_v = A_v/I*k/th;
u_v = u_v/ut;


plot(A_v(1:end-1),diff(u_v)./diff(A_v),'Color',Col,'LineWidth',3), hold on
% xlabel('$A_n$','interpreter','latex')
% ylabel('$S_n$','interpreter','latex')
xlim(xl)