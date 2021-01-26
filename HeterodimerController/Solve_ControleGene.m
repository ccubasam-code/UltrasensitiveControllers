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


% Kinetic rate of biocontroller plant
as = 20.*10^(-4)*h*.5*2;
phs = log(2)*60/3;
trs = .23*60/4;
ds = log(2)*60/30;


p2 = [as phs trs ds]; % 13-16

var = {'\kappa_c','\theta_c','\gamma_c','\delta_c','\kappa_A^+','\kappa_A^-',...
    '\kappa_R^+','\kappa_R^-','n_c','u^{tot}','\kappa_A','\kappa_I','\alpha_s',...
    '\phi_s','\psi_s','\delta_s'};
set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultAxesFontSize', 12)

p = [p1 p2 r1];

tf= 30;
tspan =[0 tf];

options = odeset('AbsTol',10^-6);
% T1 A1 R1 T2 A2 R2 T5 T6 R5 R6
x0 = [0. 0 0 0. 0 0.];

[t1,s1] = ode23s(@(t,x) ODE_ControllerGene(t,x,p),tspan,x0,options);

%
hFig=figure(10);
set(hFig,'Units','inches', 'Position', [0 100 6 3]*1.3)
%
subplot(2,1,1)
plot([0 tf],[R1 R1],'k','LineWidth',2), hold on
plot(t1,s1(:,6),'-','Color',[1 1 1]*0.8,'LineWidth',4)
hold off
subplot(2,1,2)
plot(t1,s1(:,3),'-','Color',[1 1 1]*0.8,'LineWidth',4),%hold on
% plot(t1,s1(:,4),'-','Color',[0.2 1 1]*0.8,'LineWidth',4), hold off
% ylim([0 .2])
hold off
% legend('n_s')
% k*A (p1) and th*I p(2)

%
NN = 100000;
y_c = linspace(0,1,NN)*1; %e_v(end) = 0.99999*et;
u_c = A_F_u_I(r1,y_c,p1);
idx = find (y_c<0);y_c(idx) = [];y_c(idx)=[];

u_s = logspace(-6,0,NN);
y_s = u_F_y(u_s,p2);

figure(2)
plot(y_c,u_c,'Color',[248 152 56]/255,'LineWidth',2), hold on
plot(y_s,u_s,'Color',[1 1 1]*0.5,'LineWidth',2)
plot(s1(:,6),s1(:,3),'Color',[113 111 178]/255,'LineWidth',3)
axis([0 1 0 .5])
hold off
xlabel('y'),ylabel('u')

%%

hFig=figure(3);
set(hFig,'Units','inches', 'Position', [0 100 7/2 2*1])
subplot(1,2,1)
plot([0 tf],[R1 R1],'k-','LineWidth',1), hold on
plot(t1,s1(:,6),'-','Color',[113 111 178]/255,'LineWidth',2)
xlabel('time(h)'),ylabel('y (\mu M)')
hold off
subplot(1,2,2)
plot(y_s,u_s,'Color',[160 158 158]/255,'LineWidth',2), hold on
plot(y_c,u_c,'Color',[248 152 56]/255,'LineWidth',2)
plot(s1(:,6),s1(:,3),'Color',[113 111 178]/255,'LineWidth',2)
axis([0 1 0 ut])
hold off
xlabel('y (\mu M)'),ylabel('u (\mu M)')
% axis([0 0.05 0 .5])

% Sensitive plots Control
% tf= 50;
tspan =[0 tf];
Scale = log10(2);
DIM1 = [0 6.5 7 3];
DIM2 = [0 0 7 3];
N = 5;
hFig=figure(5);
set(hFig,'Units','inches', 'Position', DIM1)
hFig=figure(6);
set(hFig,'Units','inches', 'Position', DIM2)
COLOR = [ 203 206 210;0 0 0]/255; % from red to yellow

% vec1 = [1:8]; % 9
vec1 = [10:12]; % 9
j = 1;
xx=length(vec1);
% xx = 3;
for i=vec1
    figure(5),subplot(2,4,j),plot([0 tf],[R1 R1],'-','Color',[236 24 72]/255,'LineWidth',1)
    title(var(i))
    hold on ,Plot_ScanParameter(p,i,Scale,N,tf,x0,COLOR)
    xlabel('time')
    ylim([0 1])
    
    figure(6),subplot(2,4,j),Plot_ScanParameterNCS(p,i,Scale,N,COLOR)
    xlabel('time')
    title(var(i))
    j = j+1;
end
hold off


%% Sensitive plots Process
% tf= 30;
tspan =[0 tf];
Scale = log10(2);
DIM1 = [0 6.5 7 3];
DIM2 = [0 0 7 3];
N = 5;
hFig=figure(7);
set(hFig,'Units','inches', 'Position', DIM1)
hFig=figure(8);
set(hFig,'Units','inches', 'Position', DIM2)
COLOR = [ 203 206 210;0 0 0]/255; % from red to yellow

vec1 = [13:16]; % 9
j = 1;
xx=length(vec1);
% xx = 3;
for i=vec1
    figure(7),subplot(2,4,j),plot([0 tf],[R1 R1],'-','Color',[236 24 72]/255,'LineWidth',1)
    title(var(i))
    hold on ,Plot_ScanParameter(p,i,Scale,N,tf,x0,COLOR)
    xlabel('time')
    ylim([0 1])
    
    figure(8),subplot(2,4,j),Plot_ScanParameterNCC(p,i,Scale,N,COLOR)
    xlabel('time')
    title(var(i))
    j = j+1;
end
hold off

