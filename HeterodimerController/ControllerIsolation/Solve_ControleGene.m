% Solve RNA clock
% close all,clc, clear all
% Bistable parametrs
h = 3600;
uM = 10^(-6);
uMh = uM*h;
scale = 1;
% Kinetic rate of biocontroller
kc = 20.*10^(-4)*h/4/4*4*4*scale;
thc = 20.*10^(-4)*h/4/4*4*4*scale;

gc = 3*10^4*uMh/1;
phc = log(2)*60/30; % 30 minutes
bc = 20.*10^(-4)*h/4/4;

kAp = 1*10^4*uMh/.5*0.1*10;
kAn = 10.*10^(-4)*h/1*0.1;
% kM  is 50 nM
kRp = kAp;
kRn = kAn;
ut = 0.5;
nc = 2;
KAc = 0.3; % 0.4
KIc = 0.3;

% put it in the table

r1 = 0.3; % 0.1/2/4
R1 = KIc/(thc/kc*(1+KAc/r1)-1); 

p1 = [kc thc gc phc kAp kAn kRp kRn nc ut KAc KIc]; % 1 - 12

%  100 - 1392.88/1514*100
%   100 - 757/1514*100
%   
%   [185 4884/30 4395/27 3256/20 2442/15 1953/12 1850/10]
%   [6160/50] % NEB 10-beta copmetent



var = {'\kappa_c','\theta_c','\gamma_c','\delta_c','\kappa_A^+','\kappa_A^-',...
    '\kappa_R^+','\kappa_R^-','n_c','u^{tot}','\kappa_A','\kappa_I','\alpha_s',...
    '\phi_s','\psi_s','\delta_s'};
set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultAxesFontSize', 12)



%
NN = 1000.00;
y_c = linspace(0,1,NN)*1; %e_v(end) = 0.99999*et;
u_c1 = A_F_u_I(y_c,r1,p1);
% u_c0 = A_F_u_I_V0(y_c,r1,p1); % No repressor
u_c0 = A_F_u_I_SR(y_c,r1,p1); % Single reail xR  cte
idx = find (y_c<0);y_c(idx) = [];y_c(idx)=[];

hFig=figure(2);
set(hFig,'Units','inches', 'Position', [0 6.5 3.5 1.5])

subplot(1,2,1)
plot(y_c,u_c0,'k','LineWidth',2), hold on
plot(y_c,u_c1,'b','LineWidth',2)
axis([0 1 0 .5])
hold off
xlabel('a'),ylabel('u_A')
% legend('u(x_A)','u(x_A,x_R)')
ylim([0 .5])
ax = gca;
ax.YTick = [0 .5];

%
y_c = linspace(0,1,NN)*1; %e_v(end) = 0.99999*et;
u_c1 = A_F_u_I(r1,y_c,p1);
% u_c0 = A_F_u_I_V0(r1,y_c,p1);
u_c0 = A_F_u_I_SR(r1,y_c,p1);
idx = find (y_c<0);y_c(idx) = [];y_c(idx)=[];


subplot(1,2,2)
plot(y_c,u_c0,'k','LineWidth',2), hold on
plot(y_c,u_c1,'b','LineWidth',2)
axis([0 1 0 .5])
hold off
xlabel('i'),ylabel('u_A')
% legend('u(x_A)','u(x_A,x_R)')
ylim([0 .5])
ax = gca;
ax.YTick = [0 .5];
