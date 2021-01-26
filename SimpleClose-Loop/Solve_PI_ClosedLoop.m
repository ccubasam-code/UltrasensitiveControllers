% Solve RNA clock
%close all,clc, clear all
% Bistable parametrs
h = 3600;
uM = 10^(-6);
uMh = uM*h;
F = 1;
% Kinetic rate of biocontroller
kc = 20.*10^(-4)*h/4*F;
thc = 20.*10^(-4)*h/4*1.5*F;
bc = 3*10^(4)*uMh*2*1;
ac = 3*10^(4)*uMh*2*1;
phc = log(2)*60/30; 
gc = 3*10^4*uMh;
ut = 0.5*1;

pP = [kc thc phc gc]; % 1-4
pI = [kc thc bc ac phc gc ut]; % 5-11

% Kinetic rate of biocontroller plant
as = 20.*10^(-4)*h*.5/1;
ks = .2/1; % /2
phs = log(2)*60/3;
trs = .23*60/4*1; %*2
ds = log(2)*60/30;
ns = 2;
yt = 0.6;
pS = [as ks phs trs ds ns]; % 12-17

r1 = 0.1/1.5*1; % 2 3 4
R1 = thc*r1/kc;

var = {'\kappa_c','\beta_c','\alpha_c','\phi_c','\gamma_c', 'u^{tot}', ...
    '\kappa_s','\beta_s','\alpha_s','\phi_s','\gamma_s','\theta_s','y^{tot}','w^{tot}'};
set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultAxesFontSize', 12)

c1 = 1;% proportional 
c2 = 1;% Integral
p = [pP pI pS c1 c2 r1];

tf=15;
tspan =[0 tf];

options = odeset('AbsTol',10^-6);
% T1 A1 R1 T2 A2 R2 T5 T6 R5 R6
x0 = [0 0 0 0 0 0 0 ];

[t1,s1] = ode23s(@(t,x) ODE_PI_ClosedLoop(t,x,p),tspan,x0,options);

%
hFig=figure(1);
set(hFig,'Units','inches', 'Position', [0 100 6 3]*1.3)
%
%subplot(2,3,1:3)
% plot(t0([1,end]),[r1 r1],'-','Color',[1 1 1]*0.8,'LineWidth',4)

h5 = plot(t1,s1(:,7),'Color',[113 111 178]/255,'LineWidth',3),hold on

plot(tspan,[R1 R1],'k--','LineWidth',1)
hold off
xlabel('time (h)','interpreter','latex')
ylabel('$y$ $(\mu M)$','interpreter','latex')
ax = gca;
ax.XTick = [0 tf ];
ax.YTick = [0 R1 yt];
xlim([0 tf]),ylim([0 yt])
% legend([h1 h5],'Stochastic','Deterministic')

figure(5)
subplot(3,1,1),plot(t1,s1(:,1),'b',t1,s1(:,2),'r','LineWidth',3)
legend('r_A','r_I')
subplot(3,1,2),plot(t1,s1(:,4),'b',t1,s1(:,5),'r',t1,s1(:,3),'k','LineWidth',3)
legend('r_A','r_I','u')
subplot(3,1,3),plot(t1,s1(:,6),'b',t1,s1(:,7),'r','LineWidth',3)
legend('x','y')

%%
tf=15;
tspan =[0 tf];
%Color1 = [236,226,240; 166,189,219; 28,144,153]/255;
% Color1 = [254,232,200; 253,187,132; 227,74,51]/255;
% Color1 = [102,194,165;252,141,98;141,160,203]/255;
% Color1 = [255,211,0;23,140,203;102,194,165]/255;
% Color1 = [255,241,153;162 209 234;163 218 222]/255;
Color1 = [255,231,77;93 175 219;93 191 197]/255;
Color1 = [255,231,77;93 175 219;177 209 111]/255;

options = odeset('AbsTol',10^-6);
% T1 A1 R1 T2 A2 R2 T5 T6 R5 R6
% x0 = [0 0 0 0 0 0 wt 0];


%
hFig=figure(12);
set(hFig,'Units','inches', 'Position', [0 100 6 3]*1.3)

plot(tspan,[R1 R1],'k--','LineWidth',1),hold on

% Proportional action only
p = [pP pI pS 1 0 r1];
[t1,s1] = ode23s(@(t,x) ODE_PI_ClosedLoop(t,x,p),tspan,x0,options);
plot(t1,s1(:,7),'Color',Color1(1,:),'LineWidth',3)

% Integral action only
p = [pP pI pS 0 1 r1];
[t1,s1] = ode23s(@(t,x) ODE_PI_ClosedLoop(t,x,p),tspan,x0,options);
plot(t1,s1(:,7),'Color',Color1(2,:),'LineWidth',3)

% Proportional-Integral action only
p = [pP pI pS 1 1 r1];
[t1,s1] = ode23s(@(t,x) ODE_PI_ClosedLoop(t,x,p),tspan,x0,options);
plot(t1,s1(:,7),'Color',Color1(3,:),'LineWidth',3)


legend('ref','P','I','PI')
hold off
xlabel('time (h)','interpreter','latex')
ylabel('$y$ $(\mu M)$','interpreter','latex')
ax = gca;
ax.XTick = [0 tf ];
ax.YTick = [0 R1 yt];
xlim([0 tf]),ylim([0 yt])

%% Tracking 
r1 = 0.2/1.5;
r2 = 0.1/1.5;
r3 = 0.3/1.5;
R1 = thc*r1/kc;
R2 = thc*r2/kc;
R3 = thc*r3/kc;

font_size = 12;
hFig=figure(2);
set(hFig,'Units','inches', 'Position', [0 100 3.5 2])
plot([0 tf tf 2*tf 2*tf 3*tf],[R1 R1 R2 R2 R3 R3],'k--','LineWidth',1),hold on

% Proportional action only
pS(1) = as;p = [pP pI pS 1 0 r1];
% Solve the system
[t1,s1] = ode23s(@(t,x) ODE_PI_ClosedLoop(t,x,p),tspan,x0,options);
x01 =s1(end,:);
p(end) = r2;
[t2,s2] = ode23s(@(t,x) ODE_PI_ClosedLoop(t,x,p),tspan,x01,options);
t0 = [t1 ;t2+ t1(end)];
s0 = [s1;s2];
x02 =s2(end,:);
p(end) = r3;
[t3,s3] = ode23s(@(t,x) ODE_PI_ClosedLoop(t,x,p),tspan,x02,options);
t0 = [t0 ;t3+ t0(end)];
s0 = [s0;s3];
% Plot 1
plot(t0,s0(:,7),'Color',Color1(1,:),'LineWidth',2)

% Integral action only
p = [pP pI pS 0 1 r1];
% Solve the system
[t1,s1] = ode23s(@(t,x) ODE_PI_ClosedLoop(t,x,p),tspan,x0,options);
x01 =s1(end,:);
p(end) = r2;
[t2,s2] = ode23s(@(t,x) ODE_PI_ClosedLoop(t,x,p),tspan,x01,options);
t0 = [t1 ;t2+ t1(end)];
s0 = [s1;s2];
x02 =s2(end,:);
p(end) = r3;
[t3,s3] = ode23s(@(t,x) ODE_PI_ClosedLoop(t,x,p),tspan,x02,options);
t0 = [t0 ;t3+ t0(end)];
s0 = [s0;s3];
% Plot 1
plot(t0,s0(:,7),'Color',Color1(2,:),'LineWidth',2)

% Proportional Integral action only
p = [pP pI pS 1 1 r1];
% Solve the system
[t1,s1] = ode23s(@(t,x) ODE_PI_ClosedLoop(t,x,p),tspan,x0,options);
x01 =s1(end,:);
p(end) = r2;
[t2,s2] = ode23s(@(t,x) ODE_PI_ClosedLoop(t,x,p),tspan,x01,options);
t0 = [t1 ;t2+ t1(end)];
s0 = [s1;s2];
x02 =s2(end,:);
p(end) = r3;
[t3,s3] = ode23s(@(t,x) ODE_PI_ClosedLoop(t,x,p),tspan,x02,options);
t0 = [t0 ;t3+ t0(end)];
s0 = [s0;s3];
% Plot 1
plot(t0,s0(:,7),'Color',Color1(3,:),'LineWidth',2)

hold off
ax = gca;
ax.XTick = [0 tf 2*tf 3*tf ];
ax.YTick = [0 R2  R1 R3 yt];
xlim([0 3*tf]),ylim([0 yt])

set(gca,'FontUnits','points',...
    'FontWeight','normal','FontSize',font_size,'FontName','Arial')
xlabel('time (h)','interpreter','latex','FontUnits','points',...
    'FontSize',font_size,'FontName','Arial')
ylabel('$y$ $(\mu M)$','interpreter','latex','FontUnits','points',...
    'FontSize',font_size,'FontName','Arial')
% legend({'ref','P','I','PI'},'interpreter','latex','FontUnits','points',...
%     'FontSize',font_size,'FontName','Arial','Location','NorthEast')
print -depsc2 Fig/TrackingComparation.eps
%% Disturbance Rejection 
r1 = 0.3/1.5; % 0.2

R1 = thc*r1/kc;


font_size = 12;
font_rate=10/font_size;
hFig=figure(3);
set(hFig,'Units','inches', 'Position', [0 100 3.5 2])
set(hFig,'PaperPositionMode','Auto','PaperUnits','inches','PaperSize',[3.5 2])
plot([0 3*tf],[R1 R1],'k--','LineWidth',1),hold on

% Proportional action only
pS(1) = 2*as;p = [pP pI pS 1 0 r1];
% Solve the system
[t1,s1] = ode23s(@(t,x) ODE_PI_ClosedLoop(t,x,p),tspan,x0,options);
x01 =s1(end,:);
pS(1) = 1*as;p = [pP pI pS 1 0 r1];
[t2,s2] = ode23s(@(t,x) ODE_PI_ClosedLoop(t,x,p),tspan,x01,options);
t0 = [t1 ;t2+ t1(end)];
s0 = [s1;s2];
x02 =s2(end,:);
pS(1) = 3*as;p = [pP pI pS 1 0 r1];
[t3,s3] = ode23s(@(t,x) ODE_PI_ClosedLoop(t,x,p),tspan,x02,options);
t0 = [t0 ;t3+ t0(end)];
s0 = [s0;s3];
% Plot 1
plot(t0,s0(:,7),'Color',Color1(1,:),'LineWidth',2)

% Integral action only
pS(1) = 2*as;p = [pP pI pS 0 1 r1];
% Solve the system
[t1,s1] = ode23s(@(t,x) ODE_PI_ClosedLoop(t,x,p),tspan,x0,options);
x01 =s1(end,:);
pS(1) = 1*as;p = [pP pI pS 0 1 r1];
[t2,s2] = ode23s(@(t,x) ODE_PI_ClosedLoop(t,x,p),tspan,x01,options);
t0 = [t1 ;t2+ t1(end)];
s0 = [s1;s2];
x02 =s2(end,:);
pS(1) = 3*as;p = [pP pI pS 0 1 r1];
[t3,s3] = ode23s(@(t,x) ODE_PI_ClosedLoop(t,x,p),tspan,x02,options);
t0 = [t0 ;t3+ t0(end)];
s0 = [s0;s3];
% Plot 1
plot(t0,s0(:,7),'Color',Color1(2,:),'LineWidth',2)

% Proportional Integral action only
pS(1) = 2*as;p = [pP pI pS 1 1 r1];
% Solve the system
[t1,s1] = ode23s(@(t,x) ODE_PI_ClosedLoop(t,x,p),tspan,x0,options);
x01 =s1(end,:);
pS(1) = 1*as;p = [pP pI pS 1 1 r1];
[t2,s2] = ode23s(@(t,x) ODE_PI_ClosedLoop(t,x,p),tspan,x01,options);
t0 = [t1 ;t2+ t1(end)];
s0 = [s1;s2];
x02 =s2(end,:);
pS(1) = 3*as;p = [pP pI pS 1 1 r1];
[t3,s3] = ode23s(@(t,x) ODE_PI_ClosedLoop(t,x,p),tspan,x02,options);
t0 = [t0 ;t3+ t0(end)];
s0 = [s0;s3];
% Plot 1
plot(t0,s0(:,7),'Color',Color1(3,:),'LineWidth',2)
hold off

ax = gca;
ax.XTick = [0 tf 2*tf 3*tf ];
ax.YTick = [0 R1  yt];
xlim([0 3*tf]),ylim([0 yt])
%
set(gca,'FontUnits','points',...
    'FontWeight','normal','FontSize',font_size,'FontName','Arial')
xlabel('time (h)','interpreter','latex','FontUnits','points',...
    'FontSize',font_size,'FontName','Arial')
ylabel('$y$ $(\mu M)$','interpreter','latex','FontUnits','points',...
    'FontSize',font_size,'FontName','Arial')

% legend({'ref','P','I','PI'},'interpreter','latex','FontUnits','points',...
%     'FontSize',font_size,'FontName','Arial','Location','NorthEast')
print -depsc2 Fig/DisturbanceComparation.eps