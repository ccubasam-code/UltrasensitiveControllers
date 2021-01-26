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

p1 = [kc thc bc ac phc gc ut]; % 7

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

var = {'\kappa_c','\beta_c','\alpha_c','\phi_c','\gamma_c', 'u^{tot}', ...
    '\kappa_s','\beta_s','\alpha_s','\phi_s','\gamma_s','\theta_s','y^{tot}','w^{tot}'};
set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultAxesFontSize', 12)

p = [p1 p2 r1];

tf=15;
tspan =[0 tf];

options = odeset('AbsTol',10^-6);
% T1 A1 R1 T2 A2 R2 T5 T6 R5 R6
x0 = [0 0 0 0 wt*0 0];

[t1,s1] = ode23s(@(t,x) ODE_Model_PIControler_ModuleI(t,x,p),tspan,x0,options);
x01 =s1(end,:);

% p = [p1 p2 r2];
p(end) = r2;
[t2,s2] = ode23s(@(t,x) ODE_Model_PIControler_ModuleI(t,x,p),tspan,x01,options);
t0 = [t1 ;t2+ t1(end)];
s0 = [s1;s2];
x02 =s2(end,:);
p(end) = r3;
[t3,s3] = ode23s(@(t,x) ODE_Model_PIControler_ModuleI(t,x,p),tspan,x02,options);
t0 = [t0 ;t3+ t0(end)];
s0 = [s0;s3];

%
hFig=figure(10);
set(hFig,'Units','inches', 'Position', [0 100 6 3]*1.3)
%
subplot(2,3,1:3)
plot([0 tf tf 2*tf 2*tf 3*tf],[R1 R1 R2 R2 R3 R3],'-','Color',[1 1 1]*0.8,'LineWidth',3)
hold on
h5 = plot(t0,s0(:,4),'Color',[113 111 178]/255,'LineWidth',3)
text(tf, R2*1.5,strcat('  = ',num2str(R2)),'FontSize',12);
text(2*tf, R3*1.5,strcat('  = ',num2str(R3)),'FontSize',12);
hold off
xlabel('time (h)','interpreter','latex')
ylabel('$y$ $(\mu M)$','interpreter','latex')
ax = gca;
ax.XTick = [0 tf  2*tf  3*tf ];
ax.YTick = [0 R1 yt];
xlim([0 tf*3]),ylim([0 yt])
% legend([h1 h5],'Stochastic','Deterministic')

%
N = 100000;
% e_v = linspace(0,1,N)*et;
u_v = logspace(-6,0,N)*ut; %e_v(end) = 0.99999*et;
y_v = logspace(-6,0,N)*yt;

y_e = A_F_u_I(u_v,[p1 0 r1]);
u_e = u_F_y(y_v,p2);
subplot(2,3,4),
plot(y_e,u_v,'Color',[160 158 158]/255,'LineWidth',4), hold on
plot(y_v,u_e,'Color',[248 152 56]/255,'LineWidth',1)
plot(s1(:,4),s1(:,1),'Color',[113 111 178]/255,'LineWidth',3)
hold off
xlabel('$y$ $(\mu M)$','interpreter','latex')
ylabel('$u$ $(\mu M)$','interpreter','latex')
axis([0 yt 0 ut])
ax = gca;
ax.XTick = [0  R1 yt];
ax.YTick = [0 round(s1(end,1)*10)/10 ut];

y_e = A_F_u_I(u_v,[p1 0 r2]);
% p2(1)=2*p2(1);
u_e = u_F_y(y_v,p2);
subplot(2,3,5),
plot(y_e,u_v,'Color',[160 158 158]/255,'LineWidth',4), hold on
plot(y_v,u_e,'Color',[248 152 56]/255,'LineWidth',4)
plot(s2(:,4),s2(:,1),'Color',[113 111 178]/255,'LineWidth',3)
hold off
xlabel('$y$ $(\mu M)$','interpreter','latex')
ylabel('$u$ $(\mu M)$','interpreter','latex')
axis([0 yt 0 ut])
ax = gca;
ax.XTick = [0  R2 yt];
ax.YTick = [0 round(s2(end,1)*10)/10 ut];

y_e = A_F_u_I(u_v,[p1 0 r3]);
% p2(1)=2*p2(1);
u_e = u_F_y(y_v,p2);
subplot(2,3,6),
plot(y_e,u_v,'Color',[160 158 158]/255,'LineWidth',4), hold on
plot(y_v,u_e,'Color',[248 152 56]/255,'LineWidth',4)
plot(s3(:,4),s3(:,1),'Color',[113 111 178]/255,'LineWidth',3)
hold off
xlabel('$y$ $(\mu M)$','interpreter','latex')
ylabel('$u$ $(\mu M)$','interpreter','latex')
axis([0 yt 0 ut])
ax = gca;
ax.XTick = [0  R3 yt];
ax.YTick = [0 round(s3(end,1)*10)/10 ut];

%%
print -depsc2 PerturbanceRejection.eps

%%
% p1 = [k b a ph rh ut]; % 6 Controller 
% p2 = [k th g psi ph g1t A1t]; % 7 Plan
if(1)
%     gr1 =0.4;
    p = [p1 p2 r1];
    tf= 15;
tspan =[0 tf];
COLOR = [255 172 13; 232 68 12]/255; % from red to yellow


% Parameter Scan
vec1 = [1:6]; % 9
vec2 = [7:14]; % 9
Scale = log(2);
N = 7;
DIM1 = [0 6.5 10 4];
DIM2 = [0 0 12 4];
hFig=figure(3);
set(hFig,'Units','inches', 'Position', DIM1)
j = 1;
for i=vec1
    subplot(2,3,j)
    plot([0 tf],[R1 R1],'-','Color',[236 24 72]/255,'LineWidth',3)
    hold on,Plot_ScanParameter(p,i,Scale,N,tf,x0,COLOR)
    title(var(i))
    j = j+1;
end
hold off

hFig=figure(4);
set(hFig,'Units','inches', 'Position', DIM2)
j = 1;
for i=vec2
    subplot(2,4,j)
    plot([0 tf],[R1 R1],'-','Color',[236 24 72]/255,'LineWidth',3)
    hold on ,Plot_ScanParameter(p,i,Scale,N,tf,x0,COLOR)
    title(var(i))
    j = j+1;
end
hold off
end

%% Parameter scan nullcines
DIM1 = [0 6.5 8 4];
DIM2 = [0 0 8 4];
N = 7;
hFig=figure(5);
set(hFig,'Units','inches', 'Position', DIM1)
COLOR = [255 172 13; 232 68 12]/255; % from red to yellow
vec1 = [7:10]; % 9

j = 1;
xx=length(vec1);
for i=vec1
    subplot(2,xx,j),plot([0 tf],[R1 R1],'-','Color',[153 153 153]/255,'LineWidth',3)
    title(var(i))
    hold on ,Plot_ScanParameter(p,i,Scale,N,tf,x0,COLOR)
    subplot(2,xx,j+xx),Plot_ScanParameterNCS(p,i,Scale,N,COLOR)
    title(var(i))
    j = j+1;
end
hold off

hFig=figure(6);
set(hFig,'Units','inches', 'Position', DIM2)

vec1 = [11:14]; % 9
j = 1;
xx=length(vec1);
for i=vec1
    subplot(2,xx,j),plot([0 tf],[R1 R1],'-','Color',[153 153 153]/255,'LineWidth',3)
    title(var(i))
    hold on ,Plot_ScanParameter(p,i,Scale,N,tf,x0,COLOR)
    subplot(2,xx,j+xx),Plot_ScanParameterNCS(p,i,Scale,N,COLOR)
    title(var(i))
    j = j+1;
end
hold off


%%

%% Parameter scan nullcines
Scale = log10(2);
DIM1 = [0 6.5 8 4];
DIM2 = [0 0 8 4];
N = 7;
hFig=figure(7);
set(hFig,'Units','inches', 'Position', DIM1)
COLOR = [153 153 153;0 0 0 ]/255; % from Black to grey
vec1 = [1:3]; % 9

j = 1;
xx=length(vec1);
for i=vec1
    subplot(2,xx,j),plot([0 tf],[R1 R1],'-','Color',[236 24 72]/255,'LineWidth',3)
    title(var(i))
    hold on ,Plot_ScanParameter(p,i,Scale,N,tf,x0,COLOR)
    subplot(2,xx,j+xx),Plot_ScanParameterNCC(p,i,Scale,N,COLOR)
    title(var(i))
    j = j+1;
end
hold off

hFig=figure(8);
set(hFig,'Units','inches', 'Position', DIM2)

vec1 = [4:6]; % 9
j = 1;
xx=length(vec1);
for i=vec1
    subplot(2,xx,j),plot([0 tf],[R1 R1],'-','Color',[236 24 72]/255,'LineWidth',3)
    title(var(i))
    hold on ,Plot_ScanParameter(p,i,Scale,N,tf,x0,COLOR)
    subplot(2,xx,j+xx),Plot_ScanParameterNCC(p,i,Scale,N,COLOR)
    title(var(i))
    j = j+1;
end
hold off
