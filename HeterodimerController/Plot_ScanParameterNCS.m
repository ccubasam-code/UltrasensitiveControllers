function [] = Plot_ScanParameterNCS(p,par,Scale,N,COLOR)

range = logspace(-Scale,Scale,N)*p(par);
R_v = linspace(COLOR(1,1),COLOR(2,1),N)';
B_v = linspace(COLOR(1,2),COLOR(2,2),N)';
G_v = linspace(COLOR(1,3),COLOR(2,3),N)';
Col_M = [R_v B_v G_v];

NN = 100000;

u_s = logspace(-6,0,NN);
y_s = u_F_y(u_s,p(13:16));

plot(y_s,u_s,'Color',[248 152 56]/255,'LineWidth',2), hold on

%x0 = [0 p(10) 0 0 p(11) 0 0 0];
for i=1:N
    p(par) = range(i);
    
    y_c = linspace(0,1,NN)*1; %e_v(end) = 0.99999*et;
    u_c = A_F_u_I(p(17),y_c,p);
    idx = find (y_c<0);y_c(idx) = [];y_c(idx)=[];

    plot(y_c,u_c,'Color',Col_M(i,:),'LineWidth',2)
    hold on
end
hold off
ylim([0 p(10)])
% xlim([0 p(14)])
xlim([0 1])
% xlabel('$y$ $(\mu M)$','interpreter','latex')
% ylabel('$u$ $(\mu M)$','interpreter','latex')

ax = gca;
% ax.XTick = [0  p(2)*p(16)/p(1) p(14)];
ax.XTick = [0   1];
ax.YTick = [0  p(10)];
