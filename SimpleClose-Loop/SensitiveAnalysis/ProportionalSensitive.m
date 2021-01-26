% Solve RNA clock
%close all,clc, clear all
% Bistable parametrs
h = 3600;
uM = 10^(-6);
uMh = uM*h;
% Kinetic rate of biocontroller
kc = 20.*10^(-4)*h/4*1;
thc = 20.*10^(-4)*h/4*1.*1;
bc = 3*10^4*uMh*2*1;
ac = 3*10^4*uMh*2*1;
phc = log(2)*60/30; 
gc = 3*10^4*uMh/1;
ut = 0.5;

pP = [kc thc phc gc]; % 1 - 4
pI = [kc thc bc ac phc gc ut]; % 5 - 11

% Kinetic rate of biocontroller plant
as = 20.*10^(-4)*h*.5;
ks = .2*2; %%
phs = log(2)*60/3;
trs = .23*60/4;
ds = log(2)*60/30;
ns = 2;
yt = 0.8;
pS = [as ks phs trs ds ns]; % 12-17

r1 = 0.3/1.; % .3
R1 = thc*r1/kc;

var = {'\kappa_c','\beta_c','\alpha_c','\phi_c','\gamma_c', 'u^{tot}', ...
    '\kappa_s','\beta_s','\alpha_s','\phi_s','\gamma_s','\theta_s','y^{tot}','w^{tot}'};
set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultAxesFontSize', 12)

c1 = 1;
c2 = 0;
p = [pP pI pS c1 c2 r1];

tf=15;
tspan =[0 tf];

options = odeset('AbsTol',10^-6);
% T1 A1 R1 T2 A2 R2 T5 T6 R5 R6
x0 = [0 0 0 0 0 0  0];

N = 4;
COLOR = [232 68 12 ; 255 172 13]/255; % yellow to orange
COLOR = [49 130 189; 222 235 247]/255; % light blue to blue
R_v = linspace(COLOR(2,1),COLOR(1,1),N)';
B_v = linspace(COLOR(2,2),COLOR(1,2),N)';
G_v = linspace(COLOR(2,3),COLOR(1,3),N)';
Color3 = [R_v B_v G_v];
%
font_size = 12;
hFig=figure(1);
set(hFig,'Units','inches', 'Position', [0 100 3.5 2])
%
%subplot(2,3,1:3)
% plot(t0([1,end]),[r1 r1],'-','Color',[1 1 1]*0.8,'LineWidth',4)
plot(tspan,[R1 R1],'k--','LineWidth',1),hold on

scale = [.1 1 5 10];

for i=1:N
    p(1) = pP(1)*scale(i);
    p(2) = pP(2)*scale(i);
    [t1,s1] = ode23s(@(t,x) ODE_PI_ClosedLoop(t,x,p),tspan,x0,options);
    plot(t1,s1(:,7),'Color',Color3(i,:),'LineWidth',2)
    p = [pP pI pS c1 c2 r1];
end
hold off

ax = gca;
ax.XTick = [0 tf];
ax.YTick = [0 R1 yt];
xlim([0 tf]),ylim([0 yt])

set(gca,'FontUnits','points',...
    'FontWeight','normal','FontSize',font_size,'FontName','Arial')
xlabel('time (h)','interpreter','latex','FontUnits','points',...
    'FontSize',font_size,'FontName','Arial')
ylabel('$y$ $(\mu M)$','interpreter','latex','FontUnits','points',...
    'FontSize',font_size,'FontName','Arial')
% legend({'ref','$0.1\kappa_c$','$\kappa_c$','$5\kappa_c$','$10\kappa_c$'},'interpreter','latex','FontUnits','points',...
%     'FontSize',font_size,'FontName','Arial','Location','NorthEast')
print -depsc2 Fig/PS_kappa.eps

%%

p = [pP pI pS c1 c2 r1];
hFig=figure(3);
set(hFig,'Units','inches', 'Position', [0 100 3.5 2])
%
%subplot(2,3,1:3)
% plot(t0([1,end]),[r1 r1],'-','Color',[1 1 1]*0.8,'LineWidth',4)
plot(tspan,[R1 R1],'k--','LineWidth',1),hold on

scale = [.1 1 5 10];

for i=1:N
    p(3) = pP(3)*scale(i);
    [t1,s1] = ode23s(@(t,x) ODE_PI_ClosedLoop(t,x,p),tspan,x0,options);
    plot(t1,s1(:,7),'Color',Color3(i,:),'LineWidth',2)
    p = [pP pI pS c1 c2 r1];
end
hold off

ax = gca;
ax.XTick = [0 tf];
ax.YTick = [0 R1 yt];
xlim([0 tf]),ylim([0 yt])

set(gca,'FontUnits','points',...
    'FontWeight','normal','FontSize',font_size,'FontName','Arial')
xlabel('time (h)','interpreter','latex','FontUnits','points',...
    'FontSize',font_size,'FontName','Arial')
ylabel('$y$ $(\mu M)$','interpreter','latex','FontUnits','points',...
    'FontSize',font_size,'FontName','Arial')
legend({'ref','$0.1\phi_c$','$\phi_c$','$5\phi_c$','$10\phi_c$'},'interpreter','latex','FontUnits','points',...
    'FontSize',font_size,'FontName','Arial','Location','NorthEast')
print -depsc2 Fig/PS_phi.eps

%%

p = [pP pI pS c1 c2 r1];
hFig=figure(4);
set(hFig,'Units','inches', 'Position', [0 100 3.5 2])
%
%subplot(2,3,1:3)
% plot(t0([1,end]),[r1 r1],'-','Color',[1 1 1]*0.8,'LineWidth',4)
plot(tspan,[R1 R1],'k--','LineWidth',1),hold on

scale = [.1 1 5 10];

for i=1:N
    p(4) = pP(4)*scale(i);
    [t1,s1] = ode23s(@(t,x) ODE_PI_ClosedLoop(t,x,p),tspan,x0,options);
    plot(t1,s1(:,7),'Color',Color3(i,:),'LineWidth',2)
    p = [pP pI pS c1 c2 r1];
end
hold off

ax = gca;
ax.XTick = [0 tf];
ax.YTick = [0 R1 yt];
xlim([0 tf]),ylim([0 yt])

set(gca,'FontUnits','points',...
    'FontWeight','normal','FontSize',font_size,'FontName','Arial')
xlabel('time (h)','interpreter','latex','FontUnits','points',...
    'FontSize',font_size,'FontName','Arial')
ylabel('$y$ $(\mu M)$','interpreter','latex','FontUnits','points',...
    'FontSize',font_size,'FontName','Arial')
% legend({'ref','$10^{-2}\gamma_c$','$10^{-1}\gamma_c$','$10^{0}\gamma_c$','$10^{1}\gamma_c$'},'interpreter','latex','FontUnits','points',...
%     'FontSize',font_size,'FontName','Arial','Location','NorthEast')
print -depsc2 Fig/PS_gamma.eps