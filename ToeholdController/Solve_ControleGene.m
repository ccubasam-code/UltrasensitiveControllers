% Solve RNA clock
%close all,clc, clear all
% Bistable parametrs
h = 3600;
uM = 10^(-6);
uMh = uM*h;
Scale =1;
% Kinetic rate of biocontroller

thc = 20.*10^(-4)*h/4/10*Scale;
bc = 3*10^6*uMh*2/1;
ac = 3*10^6*uMh*2/1;
phc = log(2)*60/3; % 30 minutes
gc = 3*10^4*uMh;
ut = 0.5;
KAc = 0.4;
KIc = 0.3;
trc = .23*60/4;
dc = log(2)*60/30; % 30 minutes

r1 = 0.08;% 0.2
R1 = KAc/(kc/thc*(1+KIc/r1)-1);
p1 = [kc thc bc ac phc gc ut KAc KIc]; % 1-9

% Kinetic rate of biocontroller plant
as = 20.*10^(-4)*h*.5;
ks = .2;
phs = log(2)*60/3;
trs = .23*60/4;
ds = log(2)*60/30;
ns = 2*1; %4

p2 = [as ks phs trs ds ns ]; % 8-15

var = {'\kappa_c','\theta_c','\beta_c','\alpha_c','\phi_c','\gamma_c', 'u^{tot}', ...
    '\kappa_s','\beta_s','\alpha_s','\phi_s','\gamma_s','\theta_s','y^{tot}','w^{tot}'};
set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultAxesFontSize', 12)

p = [p1 p2 r1];

tf=20;
tspan =[0 tf];

options = odeset('AbsTol',10^-6);
% T1 A1 R1 T2 A2 R2 T5 T6 R5 R6
x0 = [0. 0 0 0. 0 0];

[t1,s1] = ode23s(@(t,x) ODE_ControllerGene(t,x,p),tspan,x0,options);
x01 =s1(end,:);


%
hFig=figure(20);
set(hFig,'Units','inches', 'Position', [0 100 6 3]*1.3)
%
subplot(2,1,1)
plot([0 tf],[R1 R1],'k','LineWidth',2), hold on
plot(t1,s1(:,6),'-','Color',[1 1 1]*0.8,'LineWidth',4)
hold off
subplot(2,1,2)
plot(t1,s1(:,1)/ut,'-','Color',[1 1 1]*0.8,'LineWidth',4)
ylim([0 1])
hold off
%
N = 100000;
% e_v = linspace(0,1,N)*et;
u_v = logspace(-6,-0.001,N)*ut; %e_v(end) = [];
% y_v = logspace(-2,0,N)*trs/ds*as/phs;
y_v = logspace(-2,0,N);

y_e = A_F_u_I(u_v,[p1 0 r1]);
% y_e = I_F_u_A(u_v,[p1 r1 0]);

u_e = u_F_y(y_v,p2);
figure(3)
plot(y_e,u_v,'Color',[160 158 158]/255,'LineWidth',2), hold on
plot(y_v,u_e,'Color',[248 152 56]/255,'LineWidth',2)
plot(s1(:,6),s1(:,1),'Color',[113 111 178]/255,'LineWidth',3)
axis([0 1 0 ut])
hold off
xlabel('y'),ylabel('u')
axis([0 0.6 0 ut])

hFig=figure(3);
set(hFig,'Units','inches', 'Position', [0 100 7/2 2*1])
subplot(1,2,1)
plot([0 tf],[R1 R1],'k-','LineWidth',1), hold on
plot(t1,s1(:,6),'-','Color',[113 111 178]/255,'LineWidth',2)
xlabel('time(h)'),ylabel('y (\mu M)')
hold off
subplot(1,2,2)
plot(y_e,u_v,'Color',[160 158 158]/255,'LineWidth',2), hold on
plot(y_v,u_e,'Color',[248 152 56]/255,'LineWidth',2)
plot(s1(:,6),s1(:,1),'Color',[113 111 178]/255,'LineWidth',2)
axis([0 1 0 ut])
hold off
xlabel('y (\mu M)'),ylabel('u (\mu M)')
axis([0 0.6 0 ut])