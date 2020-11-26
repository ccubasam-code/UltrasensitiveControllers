clear all,clc
%///////////////////////////////

Omega = 600;

h = 3600;
uM = 10^(-6);
uMh = uM*h;
% Kinetic rate of biocontroller
kc = 21.*10^(-4)*h/4;
bc = 3*10^4*uMh/Omega;
ac = 3*10^4*uMh/Omega;
phc = log(2)*60/30; 
gc = 3*10^4*uMh/Omega;
ut = 0.5*Omega;

p1 = [kc bc ac phc gc ut]; % 6

% Kinetic rate of biocontroller plant
ks = kc;
bs = bc;
as = ac;
phs = phc;
gs = gc;
ths = 3.*10^(-4)*h;
yt = 0.6*Omega;
wt = yt;
r1 = 0.1*Omega;
param = 7;
p2 = [ks bs as phs gs ths yt wt];

var = {'\kappa_c','\beta_c','\alpha_c','\phi_c','\gamma_c', 'u^{tot}', ...
    '\kappa_s','\beta_s','\alpha_s','\phi_s','\gamma_s','\theta_s','y^{tot}','w^{tot}'};

set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultAxesFontSize', 12)

p = [p1 p2 r1];


op = 1;
% Volume


% Get reactions
[ S, h, endSim ] = Model_BioController( Omega,p );

% Initial condition
x0 = round(Omega*[0 0 0 0 0 0]');
% x0 = (Omega*[0.1 0 0 0 0 0  0 0]');
% Simulation time
Tmax = 15;
%///////////////////////////////


% Suggested time step
deltaT = 0.01;

% Initialize simulation:
x0 = x0;
t0 = 0;

[X1,T1,xf1,tf1] = Run_Stochastic_Sim(Omega,p,x0,t0,Tmax);

p(param) = 2*ks;
[X2,T2,xf2,tf2] = Run_Stochastic_Sim(Omega,p,xf1,tf1,2*Tmax);

p(param) = 4*ks;
[X3,T3,xf3,tf3] = Run_Stochastic_Sim(Omega,p,xf2,tf2,3*Tmax);
% Initialize "recording" variables to save the state:
% (It's much faster to reserve the whole memory at this point
% instead of enlarging the array in every step)


 X1 = X1/Omega;
 X2 = X2/Omega;
 X3 = X3/Omega;
% Plot the trajectory accumulated in arrays T and X:
% r1(1) - r4(2) - g2(3) - cr1(4) - r2(5) - r3(6) - g1(7) - cr3(8)
%
DIM1 = [0 6.5 6 3];

hFig=figure(9);
set(hFig,'Units','inches', 'Position', DIM1)

set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultAxesFontSize', 12)

% subplot(2,2,1)
% plot(T,X(1,:),'Color',[255 59 48]/255,'LineWidth',3),hold on
% plot(T,X(2,:),'Color',[0 122 255]/255,'LineWidth',3)
% plot(T,X(3,:),'Color',[255 110 30]/255,'LineWidth',3)
% % plot(T,X(4,:),'Color',[142 142 147]/255,'LineWidth',3), hold off
% legend('u','r_y','r_r')
% xlabel('time (h)')
% ylabel('\mu M')
% xlim([0 Tmax]),hold off
% ylim([0 Tmax]),ylim([0 ut/Omega])
% ax = gca;
% ax.XTick = [0 Tmax];
% ax.YTick = [0  ut/Omega];
% 
% subplot(2,2,3)
% plot(T,X(4,:),'Color',[255 59 48]/255,'LineWidth',3),hold on
% plot(T,X(5,:),'Color',[0 122 255]/255,'LineWidth',3)
% plot(T,X(6,:),'Color',[255 110 30]/255,'LineWidth',3)
% % plot(T,X(8,:),'Color',[142 142 147]/255,'LineWidth',3), hold off
% legend('y','w','z')
% xlabel('time (h)')
% ylabel('\mu M')
% hold off
% xlim([0 Tmax]),ylim([0 ut/Omega])
% ax = gca;
% ax.XTick = [0 Tmax];
% ax.YTick = [0  ut/Omega];

%
% figure(2)
% subplot(2,2,[2 4])
h1 = plot([0 3*Tmax],[r1 r1]/Omega,'-','Color',[113 111 178]/255,'LineWidth',3)
hold on
h2 = plot(T1,X1(1,:),'Color',[255 110 30]/255,'LineWidth',3);
h3 = plot(T1,X1(4,:),'Color',[1 1 1]*0.8,'LineWidth',3);
plot(T2,X2(1,:),'Color',[255 110 30]/255,'LineWidth',3);
plot(T2,X2(4,:),'Color',[1 1 1]*0.8,'LineWidth',3);
plot(T3,X3(1,:),'Color',[255 110 30]/255,'LineWidth',3);
plot(T3,X3(4,:),'Color',[1 1 1]*0.8,'LineWidth',3);
legend([h2 h3],'u','y')
xlabel('time (h)')
ylabel(' \mu M')
hold off
xlim([0 3*Tmax]),ylim([0 ut/Omega])
text(14.5, 0.4,strcat('2',var(7)),'FontSize',12);
text(29.5, 0.25,strcat('4',var(7)),'FontSize',12);
ax = gca;
ax.XTick = [0 Tmax 2*Tmax 3*Tmax];
ax.YTick = [0  r1/Omega ut/Omega];