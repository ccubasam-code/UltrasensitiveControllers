function [] = Plot_ScanParameterNCC_Break(p,par,Scale,N,COLOR)

range = logspace(-Scale,Scale,N)*p(par);
R_v = linspace(COLOR(1,1),COLOR(2,1),N)';
B_v = linspace(COLOR(1,2),COLOR(2,2),N)';
G_v = linspace(COLOR(1,3),COLOR(2,3),N)';
Col_M = [R_v B_v G_v];

NN = 2000;

for i=N:-1:1
% for i=1:N
    p(par) = range(i);
    ut = p(7)/p(8);
    u_v = logspace(-6,0,NN)*ut; %e_v(end) = 0.99999*et;
    y_e = A_F_u_I_Break(u_v,p);
    % Non-normalized plot
%     plot(y_e,u_v,'Color',Col_M(i,:),'LineWidth',4) 
    % Normalize plot
    plot(y_e*p(1)/(p(2)*p(10)),u_v/ut,'Color',Col_M(i,:),'LineWidth',4) 
%     plot(y_e,u_v/ut,'Color',Col_M(i,:),'LineWidth',4) 
    
    hold on
    
end
hold off
ylim([0 ut/ut])
xlim([0 2])
% xlabel('$y$ $(\mu M)$','interpreter','latex')
% ylabel('$u$ $(\mu M)$','interpreter','latex')

ax = gca;
% ax.XTick = [0  p(end-1) 1];
ax.YTick = [0  ut]/ut;

end
