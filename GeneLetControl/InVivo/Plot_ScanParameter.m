function [] = Plot_ScanParameter(p,par,Scale,N,tf,x0,COLOR)
range = logspace(-Scale,Scale,N)*p(par);
R_v = linspace(COLOR(1,1),COLOR(2,1),N)';
B_v = linspace(COLOR(1,2),COLOR(2,2),N)';
G_v = linspace(COLOR(1,3),COLOR(2,3),N)';
Col_M = [R_v B_v G_v];
tspan =[0 tf];
options = odeset('AbsTol',10^-6);
%x0 = [0 p(10) 0 0 p(11) 0 0 0];
for i=1:N
    p(par) = range(i);
    [t1,s1] = ode23s(@(t,x) ODE_Model_PIControler_ModuleI(t,x,p),tspan,x0,options);
    plot(t1,s1(:,4),'Color',Col_M(i,:),'LineWidth',2)
    hold on
    
end
hold off
% ylim([0 p(13)])
ylim([0 0.6])
xlim([0 tf])
xlabel('time (h)')
ylabel('$y$ $(\mu M)$','interpreter','latex')

ax = gca;
ax.XTick = [0   tf];
% ax.YTick = [0  p(15) p(13)];
ax.YTick = [0  p(15) 0.6];