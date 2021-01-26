function [] = Plot_ScanParameterNCS(p,par,Scale,N,COLOR)

range = logspace(-Scale,Scale,N)*p(par);
R_v = linspace(COLOR(1,1),COLOR(2,1),N)';
B_v = linspace(COLOR(1,2),COLOR(2,2),N)';
G_v = linspace(COLOR(1,3),COLOR(2,3),N)';
Col_M = [R_v B_v G_v];

NN = 100000;

u_v = logspace(-6,0,NN)*p(7); %e_v(end) = 0.99999*et;
y_e = A_F_u_I(u_v,[p(1:7) 0 p(16)]);
plot(y_e,u_v,'Color',[248 152 56]/255,'LineWidth',2), hold on

%x0 = [0 p(10) 0 0 p(11) 0 0 0];
for i=1:N
    p(par) = range(i);
    y_v = logspace(-6,0,NN)*p(14);
    u_e = u_F_y(y_v,p(8:15));
    plot(y_v,u_e,'Color',Col_M(i,:),'LineWidth',2)
    hold on
end
hold off
ylim([0 p(7)])
% xlim([0 p(14)])
xlim([0 0.8])
% xlabel('$y$ $(\mu M)$','interpreter','latex')
% ylabel('$u$ $(\mu M)$','interpreter','latex')

ax = gca;
% ax.XTick = [0  p(2)*p(16)/p(1) p(14)];
ax.XTick = [0   p(14)];
ax.YTick = [0  p(7)];
