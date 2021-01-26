h = 3600;
uM = 10^(-6);
uMh = uM*h;
set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultAxesFontSize', 10)
% Kinetic rate of biocontroller
kc = 20.*10^(-4)*h/4*1;
thc = 20.*10^(-4)*h/4*1;

% Update figure 3 and S4 for Single rail controller with nominal values

bc = 3*10^4*uMh*2*1;
ac = 3*10^4*uMh*2*1;
phc = log(2)*60/30; 
gc = 3*10^4*uMh;

dc = log(2)*60/30;
ut = 0.5;
rh = ut*dc;

A = 0;
I = 0.4/1; % 0.1 
w = 0.1/4;

% Fig S4, has different values of I, W

NN = 2000;
op = 1; % 0 No dilution (u cte) ---- 1 Dilution (prod/deg u)
M = 2;

p = [kc thc bc ac phc gc rh dc A I op w]; % 12
var = {'\kappa','\theta','\beta','\alpha','\phi','\gamma','\rho','\delta','a','i','op','w'};
% u_v = logspace(-6,0,NN)*ut; 
% 
% A_v1 = A_F_u_I_Break(u_v,p)*p(1)/(p(2)*p(10));
% A_v2 = A_F_u_I(u_v,p)*p(1)/(p(2)*p(10));
% 
% op = 1;
% p = [kc thc bc ac phc gc rh dc A I op w]; % 12
% A_v3 = A_F_u_I_Break(u_v,p)*p(1)/(p(2)*p(10));
% A_v4 = A_F_u_I(u_v,p)*p(1)/(p(2)*p(10));
% 
% u_v = u_v/ut;
% figure(1)
% subplot(1,2,1)
% plot(A_v1,u_v,'k','LineWidth',2), hold on
% plot(A_v2,u_v,'b','LineWidth',2), hold off
% xlim([0 1]*M)
% title('\delta_c = 0, u^{tot} = 0.5 \mu M')
% xlabel('a'), ylabel('u')
% subplot(1,2,2)
% plot(A_v3,u_v,'k','LineWidth',2), hold on
% plot(A_v4,u_v,'b','LineWidth',2), % hold off
% plot(A_v2,u_v,'b:','LineWidth',2), hold off
% xlim([0 1]*M)
% title('\delta_c > 0, u^{tot} = 0.5 \mu M')
% xlabel('a'), ylabel('u')
%% Compute sensitivity of a single parameteer
par = 4;
N = 3;
Scale = 1;
COLOR = [255 172 13; 232 68 12]/255; % from red to yellow
COLOR = [198 197 224; 113 111 178]/255; % from light purple to dark purple

figure(2)
Plot_ScanParameterNCC(p,par,Scale,N,COLOR)

%% Compute senstivity of all parameters

% op = 1;
p = [kc thc bc ac phc gc rh dc A I op w]; % 10

hFig=figure(3+op);
set(hFig,'Units','inches', 'Position', [0 6.5 7/2 3.5])

for par = 1:6
    subplot(3,3,par), Plot_ScanParameterNCC_Break(p,par,Scale,N,COLOR)
%     subplot(3,3,par), Plot_ScanParameterNCC(p,par,Scale,N,COLOR)
    title(var{par})
%     xlabel('a')
%     ylabel('u_n')
% xlim([0 2])
end
%
par = 10;
    subplot(3,3,7), Plot_ScanParameterNCC_Break(p,par,Scale,N,COLOR)
    title(var{par})
    
par = 12;
    subplot(3,3,8), Plot_ScanParameterNCC_Break(p,par,Scale,N,COLOR)
    title(var{par})

par = 7;
    subplot(3,3,9), Plot_ScanParameterNCC_Break(p,par,Scale,N,COLOR)
    title(var{par})