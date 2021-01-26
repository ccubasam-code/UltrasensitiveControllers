function [] = Compute_SensitivityInhibition(u_v,p,xl,Col)
k = p(1);
th = p(2);
ut = p(7);
u_v = u_v*ut;
A = p(8);
% Return values of titration of A, I and its hill function coefficient.
I_v = I_F_u_A(u_v,p);
I_v = I_v/A*th/k;
u_v = u_v/ut;

% Computing the hill coeficient for the activator

plot(I_v(1:end-1),-diff(u_v)./diff(I_v),'Color',Col,'LineWidth',3), hold on
% xlabel('$I_n$','interpreter','latex')
% ylabel('$n_I$','interpreter','latex')
xlim(xl)