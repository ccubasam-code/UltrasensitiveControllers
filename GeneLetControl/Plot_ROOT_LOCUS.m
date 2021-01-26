% Solve RNA clock
%close all,clc, clear all
% Bistable parametrs
h = 3600;
uM = 10^(-6);
uMh = uM*h;
scale = 1;
% Kinetic rate of biocontroller
kc = 20.*10^(-4)*h/4*scale;
thc = 20.*10^(-4)*h/4*1.5*scale;
bc = 3*10^4*uMh*2*3/3; %4
ac = 3*10^4*uMh*2*3/3; %4
phc = log(2)*60/30; 
gc = 3*10^4*uMh;
ut = 0.5;

tt1 = 2;
tt2 = 4*0;

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
r1 = 0.2;
R1 = thc*r1/kc;


p2 = [ks bs as phs gs ths yt wt];

var = {'\kappa_c','\beta_c','\alpha_c','\phi_c','\gamma_c', 'u^{tot}', ...
    '\kappa_s','\beta_s','\alpha_s','\phi_s','\gamma_s','\theta_s','y^{tot}','w^{tot}'};
set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultAxesFontSize', 12)

p = [p1 p2 r1];

tf=20;
tspan =[0 tf];

options = odeset('AbsTol',10^-6);
% T1 A1 R1 T2 A2 R2 T5 T6 R5 R6
x0 = [0 0 0 0 wt*0 0];

[t1,s1] = ode23s(@(t,x) ODE_Model_PIControler_ModuleI(t,x,p),tspan,x0,options);
x01 =s1(end,:);


% legend([h1 h5],'Stochastic','Deterministic')

%
N = 1000;
% e_v = linspace(0,1,N)*et;
u_v = logspace(-6,0,N)*ut; %e_v(end) = 0.99999*et;
y_v = logspace(-6,0,N)*yt;

y_e = A_F_u_I(u_v,[p1 0 r1]);
u_e = u_F_y(y_v,p2);



hFig=figure(1);
set(hFig,'Units','inches', 'Position', [0 100 7/2*3 2*2])
subplot(tt1,4,1+tt2)
plot([0 tf],[R1 R1],'k-','LineWidth',1), hold on
plot(t1,s1(:,4),'-','Color',[113 111 178]/255,'LineWidth',2)
xlabel('time(h)'),ylabel('y (\mu M)')
hold off
subplot(tt1,4,2+tt2)
plot(y_e,u_v,'Color',[160 158 158]/255,'LineWidth',2), hold on
plot(y_v,u_e,'Color',[248 152 56]/255,'LineWidth',2)
plot(s1(:,4),s1(:,1),'Color',[113 111 178]/255,'LineWidth',2)
axis([0 1 0 ut])
hold off
xlabel('y (\mu M)'),ylabel('u (\mu M)')
axis([0 0.6 0 ut])

%%
[yi ui] = intersections(y_e,u_v,y_v,u_e,1);
% figure(4)
% plot(y_e,u_v,'Color',[160 158 158]/255,'LineWidth',2), hold on
% plot(y_v,u_e,'Color',[248 152 56]/255,'LineWidth',2)
% plot(s1(:,4),s1(:,1),'Color',[113 111 178]/255,'LineWidth',2)
% plot(yi,ui,'ko')
% axis([0 1 0 ut])
% hold off
% xlabel('y (\mu M)'),ylabel('u (\mu M)')
% axis([0 0.6 0 ut])
% x = roots([1 2 3 4]);
%%
[rA, rI] = Compute_EQ_Contr(ui,[p1 yi r1]);
[w z] = Compute_EQ_Proc(yi,p2);
x = [rA rI ui yi w z];

% k1 = -700; k2 =200;
k1 = -220; k2 =20;
subplot(tt1,4,3+tt2)
Func_RootLocusClosedLoop(x,[p1 yi r1],p2)
% xlabel('Real Axis'),ylabel('Imag Axis'),title(' ')
axIm = findall(gcf,'String','Imaginary Axis (seconds^{-1})');
axRe = findall(gcf,'String','Real Axis (seconds^{-1})');
Bla = findall(gcf,'String','Root Locus');
set(axIm,'String','');
set(axRe,'String','');
set(Bla,'String','');
xlim([k1 k2])
ylim([-1 1]*50)

k1 = -2; k2 =1;
subplot(tt1,4,4+tt2)
Func_RootLocusClosedLoop(x,[p1 yi r1],p2)
% xlabel(' '),ylabel(' '),title(' ')
xlim([k1 k2])
ylim([-1 1]*10)
axIm = findall(gcf,'String','Imaginary Axis (seconds^{-1})');
axRe = findall(gcf,'String','Real Axis (seconds^{-1})');
Bla = findall(gcf,'String','Root Locus');
set(axIm,'String','');
set(axRe,'String','');
set(Bla,'String','');

%%
% figure(10)
% subplot(2,3,1)
% plot([0 tf],[1 1]*ui,'b--'), hold on
% plot(t1,s1(:,1),'k'), hold off
% subplot(2,3,2)
% plot([0 tf],[1 1]*rA,'b--'), hold on
% plot(t1,s1(:,2),'k'), hold off
% subplot(2,3,3)
% plot([0 tf],[1 1]*rI,'b--'), hold on
% plot(t1,s1(:,3),'k'), hold off
% %================
% subplot(2,3,4)
% plot([0 tf],[1 1]*yi,'b--'), hold on
% plot(t1,s1(:,4),'k'), hold off
% subplot(2,3,5)
% plot([0 tf],[1 1]*w,'b--'), hold on
% plot(t1,s1(:,5),'k'), hold off
% subplot(2,3,6)
% plot([0 tf],[1 1]*z,'b--'), hold on
% plot(t1,s1(:,6),'k'), hold off