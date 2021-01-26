clear all,clc
%///////////////////////////////

Omega = 1000;

h = 3600;
uM = 10^(-6);
uMh = uM*h;
% Kinetic rate of biocontroller
kc = 20.*10^(-4)*h/4;
thc = 20.*10^(-4)*h/4*1.5;
bc = 3*10^4*uMh/Omega*2;
ac = 3*10^4*uMh/Omega*2;
phc = log(2)*60/3; % 3 minutes
gc = 3*10^4*uMh/Omega;
ut = 0.5*Omega;

p1 = [kc thc bc ac phc gc ut]; % 7

% Kinetic rate of biocontroller plant
ks = kc;
bs = bc;
as = ac;
phs = phc;
gs = gc;
ths = 3.*10^(-4)*h;
yt = 0.8*Omega;
wt = yt;
r1 = 0.1*Omega;

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

p(8) = 2*ks;
[X2,T2,xf2,tf2] = Run_Stochastic_Sim(Omega,p,xf1,tf1,2*Tmax);

p(8) = 4*ks;
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
% DIM1 = [0 6.5 6 3];
DIM1 = [0 100 6 3]*1.3;
 
%
hFig=figure(10);
set(hFig,'Units','inches', 'Position', DIM1)

set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultAxesFontSize', 12)
subplot(2,3,1:3)
%
% Colort = [113 111 178]/255;
Colort = [198 197 224]/255;

% plot([0 3*Tmax],[r1 r1]/Omega,'-','Color',[102 103 102]/255,'LineWidth',4)
plot([0 3*Tmax],[r1 r1]*thc/kc/Omega,'-','Color',[102 103 102]/255,'LineWidth',4)

hold on
h1 = plot(T1,X1(4,:),'Color',Colort,'LineWidth',3);
h2 = plot(T2,X2(4,:),'Color',Colort,'LineWidth',3);
h3 = plot(T3,X3(4,:),'Color',Colort,'LineWidth',3);
% h1.Color(4) = 0.3;
% h2.Color(4) = 0.3;
% h3.Color(4) = 0.3;
