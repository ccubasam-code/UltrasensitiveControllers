function [] = Plot_ScanParameterNCS(p,par,Scale,N,COLOR)
p1 = p(1:9);
p2 = p(10:15);
r = p(16); % Reference

range = logspace(-Scale,Scale,N)*p(par);
R_v = linspace(COLOR(1,1),COLOR(2,1),N)';
B_v = linspace(COLOR(1,2),COLOR(2,2),N)';
G_v = linspace(COLOR(1,3),COLOR(2,3),N)';
Col_M = [R_v B_v G_v];

NN = 100000;
u_v = logspace(-6,-0.0003,NN)*p1(7); %e_v(end) = 0.99999*et;
y_e = A_F_u_I(u_v,[p1 0 r]);
plot(y_e,u_v,'Color',[248 152 56]/255,'LineWidth',2), hold on

%x0 = [0 p(10) 0 0 p(11) 0 0 0];
for i=1:N
    p(par) = range(i);
    p2 = p(10:15);
%     y_v = logspace(-6,0,NN)*p2(4)/p2(5)*p2(1)/p2(3)*p2(4)/p2(5);
    y_v = logspace(-6,0,NN);
    u_e = u_F_y(y_v,p2);
    idx = find ((imag(u_e)~=0)==1);
    if ~isempty(idx)
        u_e = u_e(1:idx(1));
        y_v = y_v(1:idx(1));
    end
    plot(y_v,u_e,'Color',Col_M(i,:),'LineWidth',2)
    hold on
end
hold off
ylim([0 0.5])
xlim([0 .6])
% xlabel('$y$ $(\mu M)$','interpreter','latex')
% ylabel('$u$ $(\mu M)$','interpreter','latex')

ax = gca;
ax.XTick = [0 0.12  0.6];
ax.YTick = [0  0.5];
