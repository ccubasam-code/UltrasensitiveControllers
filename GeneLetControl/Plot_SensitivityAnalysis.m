% Solve RNA clock
%close all,clc, clear all
% Bistable parametrs
h = 3600;
uM = 10^(-6);
uMh = uM*h;
% Kinetic rate of biocontroller
kc = 20.*10^(-4)*h/4;
thc = 20.*10^(-4)*h/4*1.5;
bc = 3*10^4*uMh*2;
ac = 3*10^4*uMh*2;
phc = log(2)*60/30; 
gc = 3*10^4*uMh;
ut = 0.5;

p1 = [kc thc bc ac phc gc ut]; % 6

% Kinetic rate of biocontroller plant
ks = kc;
bs = bc;
as = ac;
phs = phc;
gs = gc;
ths = 3.*10^(-4)*h;
yt = 0.8;
wt = yt;
r1 = 0.1;
r2 = 0.2;
r3 = 0.3;
R1 = thc*r1/kc;
R2 = thc*r2/kc;
R3 = thc*r3/kc;

p2 = [ks bs as phs gs ths yt wt];

var = {'\kappa_c','\theta_c','\beta_c','\alpha_c','\phi_c','\gamma_c','u^{tot}',...
    '\kappa_s','\beta_s','\alpha_s','\phi_s','\gamma_s','\theta_s','y^{tot}','w^{tot}'};

set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultAxesFontSize', 12)

p = [p1 p2 r1];

tf=10;
tspan =[0 tf];

options = odeset('AbsTol',10^-6);
% T1 A1 R1 T2 A2 R2 T5 T6 R5 R6
x0 = [0 0 0 0 wt*1 0];

%
Scale = log10(2);
DIM1 = [0 6.5 7 3];
DIM2 = [0 0 7 3];
N = 5;
hFig=figure(5);
set(hFig,'Units','inches', 'Position', DIM1)
% COLOR = [255 172 13; 232 68 12]/255; % from red to yellow
% COLOR = [ 203 206 210;77 77 77]/255; % from red to yellow
COLOR = [ 203 206 210;0 0 0]/255; % from red to yellow
% COLOR = [ 230 233 237;5 32 73]/255; % from red to yellow
vec1 = [8:11]; % 9
vec2 = [12:15]; % 9

j = 1;
xx=length(vec1);
% xx = 3;
for i=vec1
    subplot(2,xx,j),plot([0 tf],[R1 R1],'-','Color',[236 24 72]/255,'LineWidth',1)
    title(var(i))
    hold on ,Plot_ScanParameter(p,i,Scale,N,tf,x0,COLOR)
    subplot(2,xx,j+xx),Plot_ScanParameterNCS(p,i,Scale,N,COLOR)
    title(var(i))
    j = j+1;
end
hold off

hFig=figure(6);
set(hFig,'Units','inches', 'Position', DIM2)

j = 1;
% xx=length(vec2);
for i=vec2
    subplot(2,xx,j),plot([0 tf],[R1 R1],'-','Color',[236 24 72]/255,'LineWidth',1)
    title(var(i))
    hold on ,Plot_ScanParameter(p,i,Scale,N,tf,x0,COLOR)
    subplot(2,xx,j+xx),Plot_ScanParameterNCS(p,i,Scale,N,COLOR)
    title(var(i))
    j = j+1;
end
hold off
