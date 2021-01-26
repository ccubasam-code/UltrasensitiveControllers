% Solve RNA clock
%close all,clc, clear all
% Bistable parametrs
h = 3600;
uM = 10^(-6);
uMh = uM*h;
scale = 10;
% Kinetic rate of biocontroller
kc = 20.*10^(-4)*h/4*scale;
thc = 20.*10^(-4)*h/4*1.5*scale;
bc = 3*10^4*uMh*2;
ac = 3*10^4*uMh*2;
phc = log(2)*60/30; 
gc = 3*10^4*uMh;
ut = 0.5;

p1 = [kc thc bc ac phc gc ut]; % 6

% Kinetic rate of biocontroller plant
ks = 20.*10^(-4)*h/4;
bs = 3*10^4*uMh*2;
as = 3*10^4*uMh*2;
phs = log(2)*60/30;
gs = 3*10^4*uMh;
ths = 3.*10^(-4)*h;

yt = 0.8;
wt = yt;
r1 = 0.2/2*4;
R1 = thc*r1/kc;


p2 = [ks bs as phs gs ths yt wt];

var = {'\kappa_c','\beta_c','\alpha_c','\phi_c','\gamma_c', 'u^{tot}', ...
    '\kappa_s','\beta_s','\alpha_s','\phi_s','\gamma_s','\theta_s','y^{tot}','w^{tot}'};
set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultAxesFontSize', 12)

p = [p1 p2 r1];

%%%
p = [1.8000    2.7000  216.0000  216.0000    1.3863  108.0000    1.0000    3.4510 ...
    199.3512  245.0301    0.8084   71.2388    0.6769    0.4210    0.5838    0.3296];
r1 = p(end);
%%%

tf=20;
tspan =[0 tf];

options = odeset('AbsTol',10^-6);
% T1 A1 R1 T2 A2 R2 T5 T6 R5 R6
x0 = [0 0 0 0 wt*0 0];

[t1,s1] = ode23s(@(t,x) ODE_Model_PIControler_ModuleI(t,x,p),tspan,x0,options);
x01 =s1(end,:);


% legend([h1 h5],'Stochastic','Deterministic')

%
NN = 1000;

u_c = logspace(-6,0,NN)*p(7); %e_v(end) = 0.99999*et;
y_c = A_F_u_I(u_c,[p(1:7) 0 r1]);

y_s = linspace(0,1,NN)*p(14);
u_s = u_F_y(y_s,p(8:15));
idx = find (y_c<0);u_c(idx) = [];y_c(idx)=[];

[u_i,y_i] = intersections(u_c,y_c,u_s,y_s,1);
        


hFig=figure(3);
set(hFig,'Units','inches', 'Position', [0 100 7/2 2*1])
subplot(1,2,1)
plot([0 tf],[R1 R1],'k-','LineWidth',1), hold on
plot(t1,s1(:,4),'-','Color',[113 111 178]/255,'LineWidth',2)
xlabel('time(h)'),ylabel('y (\mu M)')
hold off
subplot(1,2,2)
plot(y_c,u_c,'Color',[248 152 56]/255,'LineWidth',2), hold on
plot(y_s,u_s,'Color',[160 158 158]/255,'LineWidth',2)
plot(y_i,u_i,'ko','LineWidth',2)
hold off
ylim([0 1]*ut)

%
[ExistEquilibrium,eigJ,eq] = ComputeStability(p,u_i,y_i);
eigJ
eq

figure(4)
ii=1;
subplot(2,3,ii), plot(t1,s1(:,ii),'k'), hold on
plot(tf,eq(ii),'ro'), hold off
ylabel('u')
ii=ii+1;
subplot(2,3,ii), plot(t1,s1(:,ii),'k'), hold on
plot(tf,eq(ii),'ro'), hold off
ylabel('r_A')
ii=ii+1;
subplot(2,3,ii), plot(t1,s1(:,ii),'k'), hold on
plot(tf,eq(ii),'ro'), hold off
ylabel('r_I')
ii=ii+1;
subplot(2,3,ii), plot(t1,s1(:,ii),'k'), hold on
plot(tf,eq(ii),'ro'), hold off
ylabel('y')
ii=ii+1;
subplot(2,3,ii), plot(t1,s1(:,ii),'k'), hold on
plot(tf,eq(ii),'ro'), hold off
ylabel('w')
ii=ii+1;
subplot(2,3,ii), plot(t1,s1(:,ii),'k'), hold on
plot(tf,eq(ii),'ro'), hold off
ylabel('z')