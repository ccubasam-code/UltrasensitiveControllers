% Solve RNA clock
%close all,clc, clear all
% Bistable parametrs
h = 3600;
uM = 10^(-6);
uMh = uM*h;
zz = 3;
F = 4;
SC = 1; % 100
% Kinetic rate of biocontroller
kc = 20.*10^(-4)*h/4*F;
thc = 20.*10^(-4)*h/4*1.*F;
bc = 3*10^(4)*uMh*2*SC;  % 1 - 10 - 100
ac = 3*10^(4)*uMh*2*SC;  % 1 - 10 - 100
phc = log(2)*60/30; % 1 1/5 1/10 
gc = 3*10^4*uMh*1;
ut = 0.5*1; % 1 5 10
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

r1 = 0.1/1.*1; % 2 3 4
R1 = thc*r1/kc;

var = {'\kappa_c','\beta_c','\alpha_c','\phi_c','\gamma_c', 'u^{tot}', ...
    '\kappa_s','\beta_s','\alpha_s','\phi_s','\gamma_s','\theta_s','y^{tot}','w^{tot}'};
set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultAxesFontSize', 12)

c1 = 0;% proportional 
c2 = 1;% Integral
p = [pP pI pS c1 c2 r1];

tf=15;
tspan =[0 tf];

options = odeset('AbsTol',10^-6);
% T1 A1 R1 T2 A2 R2 T5 T6 R5 R6
x0 = [0 0 0 0 0 0 0 ];

Color1 = [255,231,77;93 175 219;177 209 111]/255;

%% Tracking 
r1 = 0.2/1.;
r2 = 0.1/1.;
r3 = 0.3/1.;

R1 = thc*r1/kc;
R2 = thc*r2/kc;
R3 = thc*r3/kc;

font_size = 12;
hFig=figure(2);
set(hFig,'Units','inches', 'Position', [0 100 7 2*3])
subplot(3,2,2*(zz-1)+1)
plot([0 tf tf 2*tf 2*tf 3*tf],[R1 R1 R2 R2 R3 R3],'k--','LineWidth',1),hold on

% Proportional action only
pS(1) = as;p = [pP pI pS 1 0 r1];
% Solve the system
[t1,s1] = ode23s(@(t,x) ODE_Integral_ClosedLoop(t,x,p),tspan,x0,options);
x01 =s1(end,:);
p(end) = r2;
[t2,s2] = ode23s(@(t,x) ODE_Integral_ClosedLoop(t,x,p),tspan,x01,options);
t0 = [t1 ;t2+ t1(end)];
s0 = [s1;s2];
x02 =s2(end,:);
p(end) = r3;
[t3,s3] = ode23s(@(t,x) ODE_Integral_ClosedLoop(t,x,p),tspan,x02,options);
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


legend('r','MS','BM')
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

%% Tracking with non-ideal controller

%% Tracking 
r1 = 0.2/1.;
r2 = 0.1/1.;
r3 = 0.3/1.;

R1 = thc*r1/kc;
R2 = thc*r2/kc;
R3 = thc*r3/kc;

font_size = 12;
% hFig=figure(3);
% set(hFig,'Units','inches', 'Position', [0 5 3.5 2])
subplot(3,2,2*(zz-1)+2)
plot([0 tf tf 2*tf 2*tf 3*tf],[R1 R1 R2 R2 R3 R3],'k--','LineWidth',1),hold on

% Proportional action only
pS(1) = as;p = [pP pI pS 1 0 r1];
% Solve the system
[t1,s1] = ode23s(@(t,x) ODE_Integral_ClosedLoopDycU(t,x,p),tspan,x0,options);
x01 =s1(end,:);
p(end) = r2;
[t2,s2] = ode23s(@(t,x) ODE_Integral_ClosedLoopDycU(t,x,p),tspan,x01,options);
t0 = [t1 ;t2+ t1(end)];
s0 = [s1;s2];
x02 =s2(end,:);
p(end) = r3;
[t3,s3] = ode23s(@(t,x) ODE_Integral_ClosedLoopDycU(t,x,p),tspan,x02,options);
t0 = [t0 ;t3+ t0(end)];
s0 = [s0;s3];
% Plot 1
plot(t0,s0(:,7),'Color',Color1(1,:),'LineWidth',2)

% Integral action only
p = [pP pI pS 0 1 r1];
% Solve the system
[t1,s1] = ode23s(@(t,x) ODE_Integral_ClosedLoopDycU(t,x,p),tspan,x0,options);
x01 =s1(end,:);
p(end) = r2;
[t2,s2] = ode23s(@(t,x) ODE_Integral_ClosedLoopDycU(t,x,p),tspan,x01,options);
t0 = [t1 ;t2+ t1(end)];
s0 = [s1;s2];
x02 =s2(end,:);
p(end) = r3;
[t3,s3] = ode23s(@(t,x) ODE_Integral_ClosedLoopDycU(t,x,p),tspan,x02,options);
t0 = [t0 ;t3+ t0(end)];
s0 = [s0;s3];
% Plot 1
plot(t0,s0(:,7),'Color',Color1(2,:),'LineWidth',2)

hold off

% legend('r','MS','BM')

ax = gca;
ax.XTick = [0 tf 2*tf 3*tf ];
ax.YTick = [0 R2  R1 R3 yt];
xlim([0 3*tf]),ylim([0 yt])

% legend('ref','MS','Brink')
set(gca,'FontUnits','points',...
    'FontWeight','normal','FontSize',font_size,'FontName','Arial')
xlabel('time (h)','interpreter','latex','FontUnits','points',...
    'FontSize',font_size,'FontName','Arial')
ylabel('$y$ $(\mu M)$','interpreter','latex','FontUnits','points',...
    'FontSize',font_size,'FontName','Arial')