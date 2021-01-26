h = 3600;
uM = 10^(-6);
uMh = uM*h;
set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultAxesFontSize', 10)
% Kinetic rate of biocontroller
kc = 20.*10^(-4)*h/4*1;
thc = 20.*10^(-4)*h/4*1;

bc = 3*10^4*uMh*2*1;
ac = 3*10^4*uMh*2*1;
phc = log(2)*60/30; 
gc = 3*10^4*uMh;

dc = log(2)*60/30;
ut = 0.5;
rh = ut*dc;

% To do plots for paper.
% Check the plots for rho and delta to make it conistant ut=0.5


A = 0;
I = 0.4;

NN = 1000;
op = 1; % 0: u constat - 1: production rates
w = 0;

p = [kc thc bc ac phc gc rh dc A I op w]; % 10
var = {'\kappa','\theta','\beta','\alpha','\phi','\gamma','\rho','\delta','a','i'};
u_v = logspace(-6,0,NN)*ut; 

A_v = A_F_u_I(u_v,p);

figure(1)

plot(A_v,u_v,'k','LineWidth',2)


%% Compute sensitivity of a single parameteer
par = 4;
N = 3;
Scale = 1;
COLOR = [255 172 13; 232 68 12]/255; % from yellow to red
% COLOR = [0 0 0; 0 0 1]; % from red to yellow

figure(2)
Plot_ScanParameterNCC(p,par,Scale,N,COLOR)

%% Compute senstivity of all parameters

hFig=figure(10+op);
set(hFig,'Units','inches', 'Position', [0 6.5 7/2 3.5])

for par = 1:8
    subplot(3,3,par), Plot_ScanParameterNCC(p,par,Scale,N,COLOR)
    title(var{par})
%     xlabel('a')
%     ylabel('u_n')
end
