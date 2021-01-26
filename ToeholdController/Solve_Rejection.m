% Solve RNA clock
%close all,clc, clear all
% Bistable parametrs
h = 3600;
uM = 10^(-6);
uMh = uM*h;
% Kinetic rate of biocontroller
kc = 20.*10^(-4)*h/4/10;
thc = 20.*10^(-4)*h/4/10;
bc = 3*10^6*uMh*2;
ac = 3*10^6*uMh*2;
phc = log(2)*60/3; % 30 minutes
gc = 3*10^4*uMh;
ut = 0.5;
KAc = 0.4;
KIc = 0.3;
trc = .23*60/4;
dc = log(2)*60/30; % 30 minutes

r1 = 0.09;
r2 = 0.18;
r3 = 0.27;

R1 = KAc/(kc/thc*(1+KIc/r1)-1); 
R2 = KAc/(kc/thc*(1+KIc/r2)-1); 
R3 = KAc/(kc/thc*(1+KIc/r3)-1); 
p1 = [kc thc bc ac phc gc ut KAc KIc]; % 1-9

% Kinetic rate of biocontroller plant
as = 20.*10^(-4)*h*.5;
ks = .2;
phs = log(2)*60/3;
trs = .23*60/4;
ds = log(2)*60/30;
ns = 2;

nn = 10;
p2 = [as ks phs trs ds ns]; % 10-15

var = {'\kappa_c','\theta_c','\beta_c','\alpha_c','\phi_c','\gamma_c', 'u^{tot}', ...
    '\kappa_s','\beta_s','\alpha_s','\phi_s','\gamma_s','\theta_s','y^{tot}','w^{tot}'};
set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultAxesFontSize', 12)

p = [p1 p2 r1];
pp = p;
tf=50;
tspan =[0 tf];

options = odeset('AbsTol',10^-6);
% T1 A1 R1 T2 A2 R2 T5 T6 R5 R6
x0 = [0. 0 0 0. 0 0];

[t1,s1] = ode23s(@(t,x) ODE_ControllerGene(t,x,p),tspan,x0,options);
x01 =s1(end,:);

p(nn) = 2*pp(nn);
[t2,s2] = ode23s(@(t,x) ODE_ControllerGene(t,x,p),tspan,x01,options);
t0 = [t1 ;t2+ t1(end)];
s0 = [s1;s2];
x02 =s2(end,:);
p(nn) = 4*pp(nn);
[t3,s3] = ode23s(@(t,x) ODE_ControllerGene(t,x,p),tspan,x02,options);
t0 = [t0 ;t3+ t0(end)];
s0 = [s0;s3];
%
hFig=figure(10);
set(hFig,'Units','inches', 'Position', [0 100 6 3]*1.3)
%
subplot(2,3,1:3)
plot([0  3*tf],[R1 R1],'k','LineWidth',2), hold on
plot(t0,s0(:,6),'-','Color',[113 111 178]/255,'LineWidth',3)

hold off
xlabel('time (h)','interpreter','latex')
ylabel('$y$ $(\mu M)$','interpreter','latex')
ax = gca;
ax.XTick = [0 tf  2*tf  3*tf ];
ax.YTick = [0 round(100*R1)/100 round(YT*10)/10];
xlim([0 tf*3]),ylim([0 round(YT*10)/10])


%
% p1 = [kc thc bc ac phc gc ut KAc KIc]; % 1-9
% p2 = [as ks phs trs ds ns]; % 10-15
% p = [p1 p2 r];
p(nn) = 1*pp(nn);p1 = p(1:9);p2 = p(10:15);

N = 10000;
u_v = logspace(-6,-0.0003,N)*ut; %e_v(end) = [];
y_v = logspace(-2,0,N)*trs/ds*as/phs;
y_e = A_F_u_I(u_v,[p1 0 r1]);
u_e = u_F_y(y_v,p2);
YT = as/phs*trs/ds;
UT = ut;
subplot(2,3,4),
plot(y_e,u_v,'Color',[248 152 56]/255,'LineWidth',4), hold on
plot(y_v,u_e,'Color',[160 158 158]/255,'LineWidth',4)
plot(s1(:,6),s1(:,1),'Color',[113 111 178]/255,'LineWidth',3)
hold off
xlabel('$y$ $(\mu M)$','interpreter','latex')
ylabel('$u$ $(\mu M)$','interpreter','latex')
axis([0 YT 0 UT])
ax = gca;

ax.XTick = [0  round(R1*100)/100 round(YT*10)/10];
ax.YTick = [0 round(s1(end,1)*10)/10 round(UT*10)/10];

%
u_v = logspace(-6,-0.0003,N)*ut; %e_v(end) = [];
y_v = logspace(-2,0,N)*trs/ds*as/phs;

p(nn) = 2*pp(nn);p1 = p(1:9);p2 = p(10:15);

y_e = A_F_u_I(u_v,[p1 0 r1]);
% p2(1) = 2*pp(nn);
u_e = u_F_y(y_v,p2);

subplot(2,3,5),
plot(y_e,u_v,'Color',[248 152 56]/255,'LineWidth',4), hold on
plot(y_v,u_e,'Color',[160 158 158]/255,'LineWidth',4)
plot(s2(:,6),s2(:,1),'Color',[113 111 178]/255,'LineWidth',3)
hold off
xlabel('$y$ $(\mu M)$','interpreter','latex')
ylabel('$u$ $(\mu M)$','interpreter','latex')
axis([0 YT 0 UT])
ax = gca;
ax.XTick = [0  round(R1*100)/100 round(YT*10)/10];
ax.YTick = [0 round(s2(end,1)*10)/10 round(UT*10)/10];

%
u_v = logspace(-6,-0.0003,N)*ut; %e_v(end) = [];
y_v = logspace(-2,0,N)*trs/ds*as/phs;

p(nn) = 4*pp(nn);p1 = p(1:9);p2 = p(10:15);
y_e = A_F_u_I(u_v,[p1 0 r1]);
% p2(1) = 4*pp(nn);
u_e = u_F_y(y_v,p2);

subplot(2,3,6),
plot(y_e,u_v,'Color',[248 152 56]/255,'LineWidth',4), hold on
plot(y_v,u_e,'Color',[160 158 158]/255,'LineWidth',4)
plot(s3(:,6),s3(:,1),'Color',[113 111 178]/255,'LineWidth',3)
hold off
xlabel('$y$ $(\mu M)$','interpreter','latex')
ylabel('$u$ $(\mu M)$','interpreter','latex')
axis([0 YT 0 UT])
ax = gca;
ax.XTick = [0  round(R1*100)/100 round(YT*10)/10];
ax.YTick = [0 round(s3(end,1)*10)/10 round(UT*10)/10];
