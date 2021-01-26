% Solve RNA clock
%close all,clc, clear all
% Bistable parametrs
h = 3600;
uM = 10^(-6);
uMh = uM*h;
set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultAxesFontSize', 10)
% Kinetic rate of biocontroller
kc = 20.*10^(-4)*h/4;
thc = 20.*10^(-4)*h/4*1.5;
bc = 3*10^4*uMh*2;
ac = 3*10^4*uMh*2;
phc = log(2)*60/30; % 3 minutes 
gc = 3*10^4*uMh;
ut = 0.5;

A = 0.1;
I = 0.1;
p = [kc thc bc ac phc gc ut A I]; % 6

var = {'\kappa_c','\theta_c','\beta_c','\alpha_c','\phi_c','\gamma_c','u^{tot}','A^{tot}','I^{tot}'};

N = 100000;
% e_v = linspace(0,1,N)*et;
u_v = logspace(-6,0,N)*ut; %e_v(end) = 0.99999*et;

I_v = I_F_u_A(u_v,p);
A_v = A_F_u_I(u_v,p);

Color3 = [8 104 172; 67 162 202; 123 204 196]/255;
% Computing derivative


% Normalization
u_v = u_v/ut;
I_v = I_v/A*thc/kc;
A_v = A_v/I*kc/thc;

zh = [0 2];
figure(1)
subplot(1,2,1),plot(I_v,u_v,'k','LineWidth',3)
xlim(zh),
xlabel('I/A'), ylabel('u/u^{tot}')
subplot(1,2,2),plot(A_v,u_v,'k','LineWidth',3)
xlim(zh),
xlabel('A/I'), ylabel('u/u^{tot}')

%%
close all
DIM1 = [0 0 3.5+0.5 2];
hFig=figure(2);
set(gcf,'PaperPositionMode','auto')
set(hFig,'Units','inches', 'Position', DIM1)
figure(2)
subplot(2,2,1)
plot(A_v,u_v,'k','LineWidth',3)
xlim(zh),
xlabel('I_n'), ylabel('u_n')

ax = gca;ax.YTick = [0 1];

subplot(2,2,3)

plot(A_v(1:end-1),diff(u_v)./diff(A_v),'r','LineWidth',3)
xlim(zh),
xlabel('A_n'), ylabel('S_n')
ylim([0 12])
% ax = gca;ax.YTick = [0 10];

subplot(2,2,2)
plot(I_v,u_v,'k','LineWidth',3)
xlim(zh),
xlabel('I_n'), ylabel('u_n')
ax = gca;ax.YTick = [0 1];
subplot(2,2,4)

plot(I_v(1:end-1),-diff(u_v)./diff(I_v),'r','LineWidth',3)
xlim(zh),
xlabel('I_n'), ylabel('S_n')
ylim([0 12])
ax = gca;ax.YTick = [0 10];

% print -depsc2 test1.eps
% set(gcf, 'PaperUnits', 'inches');
% set(gcf, 'PaperPosition', DIM1); %
% set(gcf, 'PaperSize', DIM1); %
% % 
% saveas(gcf,['test2.eps'],'epsc')
%%

scale3 = [10 .1  0.01];
DIM1 = [0 60 7.25 2];

hFig=figure(2);
set(hFig,'Units','inches', 'Position', DIM1)
% A_v = A_v/I*kc/thc;
j=3;
N = 100000;
% e_v = linspace(0,1,N)*et;
u_v = logspace(-6,0,N); %e_v(end) = 0.99999*et;
op = 0; % only for 8
xx = 1;
xlm = [0 2.2]
for z=1:3
        p0 = p;
        p0(j) = p(j)*scale3(z);
%         subplot(2,4,xx),
        subplot(1,2,xx),
        if (op==1)
            I_v = I_F_u_A(u_v*p0(7),p0);
%             I_v = I_v/p0(8)*p0(2)/p0(1);
            plot(I_v,u_v,'Color',Color3(z,:),'LineWidth',3), hold on
        else
            A_v = A_F_u_I(u_v*p0(7),p0);
%             A_v = A_v/p0(9)*p0(1)/p0(2);
            plot(A_v,u_v,'Color',Color3(z,:),'LineWidth',3), hold on
        end
end
Func_AddLegend(0.22,0.5,j,0.02,0.2,scale3,var,Color3)
% Func_AddLegend(.65,0.6,j,0.05,0.25,scale3,var,Color3)
xlim([0 0.3])
ax = gca;
% ax.XTick = [ 0.1   0.3 0.6];
ax.XTick = [0 0.15 0.3];
ax.YTick = [0 1];
hold off
for z=1:3
        p0 = p;
        p0(j) = p(j)*scale3(z);
        if op ==1
%             subplot(2,4,xx+1),Compute_SensitivityInhibition(u_v,p0, xlm,Color3(z,:))
            subplot(1,2,xx+1),Compute_SensitivityInhibition(u_v,p0, xlm,Color3(z,:))
        else
%             subplot(2,4,xx+1),Compute_SensitivityActivation(u_v,p0, xlm,Color3(z,:))
            subplot(1,2,xx+1),Compute_SensitivityActivation(u_v,p0, xlm,Color3(z,:))
        end
        ylim([0 yl(j)])
end
Func_AddLegend(x(j),y(j),j,aa(j),bb(j),scale3,var,Color3)
ax = gca;
ax.XTick = [0 1 2];
% ax.YTick = [0 10];
%%
scale3 = [8 4  1.5];
DIM1 = [0 60 7.25 2];

hFig=figure(2);
set(hFig,'Units','inches', 'Position', DIM1)
% A_v = A_v/I*kc/thc;
j=8;
N = 100000;
% e_v = linspace(0,1,N)*et;
u_v = logspace(-6,0,N); %e_v(end) = 0.99999*et;
op = 1; % only for 8
xx = 3;
for z=1:3
        p0 = p;
        p0(j) = p(j)*scale3(z);
        subplot(2,4,xx),
        if (op==1)
            I_v = I_F_u_A(u_v*p0(7),p0);
%             I_v = I_v/p0(8)*p0(2)/p0(1);
            plot(I_v,u_v,'Color',Color3(z,:),'LineWidth',3), hold on
        else
            A_v = A_F_u_I(u_v*p0(7),p0);
%             A_v = A_v/p0(9)*p0(1)/p0(2);
            plot(A_v,u_v,'Color',Color3(z,:),'LineWidth',3), hold on
        end
end
% Func_AddLegend(0.22,0.5,j,0.02,0.2,scale3,var,Color3)
Func_AddLegend(.65,0.6,j,0.05,0.25,scale3,var,Color3)
xlim([0 1.05])
ax = gca;
ax.XTick = [ 0.1   0.3 0.6];
% ax.XTick = [0 0.15 0.3];
ax.YTick = [0 1];
hold off
for z=1:3
        p0 = p;
        p0(j) = p(j)*scale3(z);
        if op ==1
            subplot(2,4,xx+1),Compute_SensitivityInhibition(u_v,p0, xlm,Color3(z,:))
        else
            subplot(2,4,xx+1),Compute_SensitivityActivation(u_v,p0, xlm,Color3(z,:))
        end
        ylim([0 yl(j)])
end
Func_AddLegend(x(j),y(j),j,aa(j),bb(j),scale3,var,Color3)
ax = gca;
% ax.XTick = [0 1 2];
ax.YTick = [0 10];
%% 


N = 100000;
% e_v = linspace(0,1,N)*et;
u_v = logspace(-6,0,N); %e_v(end) = 0.99999*et;


scale3 = [10  1  .1];
% Color3 = [232 68 12;255 102 0; 255 172 13]/255; % red orange yellow
% Color3 = [0 109 44;44 162 95; 102 194 164]/255; % scale green
Color3 = [8 104 172; 67 162 202; 123 204 196]/255;

x = [1.5 1.5 1.5 1.5 1.45 1.45 1.35 1.35 1.4];
y = [14 14 36 36 64*2 14 62*2 14 14]/2;
aa = [0.15 0.15 0.15 0.15 0.15 0.15 0.15 0.15 0.15];

xlm = [0 2.2]
DIM1 = [0 0 7.25 2];

hFig=figure(1);
set(hFig,'Units','inches', 'Position', DIM1)

yl = [12 12 30 30 110 12 110 12 12];
tt = [10 10 30 30 100 10 100 10 10];
bb = y/2.5;
for j= [1:7]
    p0 = p;
    if( j~=5)
        for z=1:3
            p0(j) = p(j)*scale3(z);
            subplot(2,4,j),Compute_SensitivityActivation(u_v,p0, xlm,Color3(z,:))   
        end
    else
        for z=3:-1:1
            p0(j) = p(j)*scale3(z);
            subplot(2,4,j),Compute_SensitivityActivation(u_v,p0, xlm,Color3(z,:))   
        end
    end
    Func_AddLegend(x(j),y(j),j,aa(j),bb(j),scale3,var,Color3)
    ylim([0 yl(j)])
    ax = gca;
    ax.YTick = [0 tt(j)];
end
j = 9;
for z=1:3
        p0 = p;
        p0(j) = p(j)*scale3(z);
        subplot(2,4,j-1),Compute_SensitivityInhibition(u_v,p0, xlm,Color3(z,:))
        ylim([0 yl(j)])
end
Func_AddLegend(x(j),y(j),j,aa(j),bb(j),scale3,var,Color3)
ax = gca;
ax.YTick = [0 tt(j)];

%%
%% 


N = 100000;
% e_v = linspace(0,1,N)*et;
u_v = logspace(-6,0,N); %e_v(end) = 0.99999*et;


scale3 = [10  1  .1];
% Color3 = [232 68 12;255 102 0; 255 172 13]/255; % red orange yellow
% Color3 = [0 109 44;44 162 95; 102 194 164]/255; % scale green
Color3 = [8 104 172; 67 162 202; 123 204 196]/255;

x = [1.5 1.5 1.5 1.5 1.45 1.45 1.35 1.35 1.4]+0.1;
y = [1 1 1 0.8 .8 1 0.8 1 1]/2.;
aa = [0.15 0.15 0.15 0.15 0.15 0.15 0.15 0.15 0.15];

xlm = [0 2.2]
DIM1 = [0 10 7.25 3];

hFig=figure(1);
set(hFig,'Units','inches', 'Position', DIM1)

yl = [1 1 1 1 1 1 1 1 1];
tt = [10 10 30 30 100 10 100 10 10];
bb = [1 1 1 1 1 1 1 1 1]/3.8;
xxlim = [2 2 0.3 0.3 0.4 0.3 0.3 2 2];
jj=1;
for j= [1:7,9]
    p0 = p;
    for z=1:3
        p0(j) = p(j)*scale3(z);
        A_v = A_F_u_I(u_v*p0(7),p0);
            A_v = A_v/p0(9)*p0(1)/p0(2);
            subplot(3,3,jj),plot(A_v,u_v,'Color',Color3(z,:),'LineWidth',3), hold on
%             xlim([0 xxlim(j)])
    end
    Func_AddLegend(x(j),y(j),j,aa(j),bb(j),scale3,var,Color3)
    ylim([0 yl(j)])
    xlim([0 2])
    ax = gca;
    ax.YTick = [0 1];
    ax.XTick = [0 1 2];
    jj = jj+1;
end
