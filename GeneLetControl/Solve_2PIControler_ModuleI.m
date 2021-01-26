% Solve RNA clock
%close all,clc, clear all
% Bistable parametrs
h = 3600;
uM = 10^(-6);
uMh = uM*h;
% Kinetic rate of biocontroller
kc = 21.*10^(-4)*h/4;
bc = 3*10^4*uMh;
ac = 3*10^4*uMh;
phc = log(2)*60/30; 
gc = 3*10^4*uMh;
ut = 0.5;

p1 = [kc bc ac phc gc ut]; % 6

p2 = [kc bc ac phc gc ut]; % 6

% Kinetic rate of biocontroller plant
ks = kc;
bs = bc;
as = ac;
phs = phc;
gs = gc;
ths = 3.*10^(-4)*h;
yt = 0.6;
wt = yt;
r1 = .01;
r2 = 0.25;
p3 = [ks bs as phs gs ths yt wt];

var = {'\kappa_c','\beta_c','\alpha_c','\phi_c','\gamma_c', 'u^{tot}', ...
    '\kappa_s','\beta_s','\alpha_s','\phi_s','\gamma_s','\theta_s','y^{tot}','w^{tot}'};
set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultAxesFontSize', 12)

p = [p1 p2 p3 r1 r2];

tf=15*1;
tspan =[0 tf];

options = odeset('AbsTol',10^-6);
% T1 A1 R1 T2 A2 R2 T5 T6 R5 R6
x0 = [0.2 0 0 0 0 0 0 wt 0];

[t1,s1] = ode23s(@(t,x) ODE_Model_2PIControler_ModuleI(t,x,p),tspan,x0,options);
x01 =s1(end,:);

% t0 = t1; s0 = s1;
p = [p1 p2 p3 r1 r2];
ii = 3;
fac = 1;
p(ii) = p(ii)*fac;
[t2,s2] = ode23s(@(t,x) ODE_Model_2PIControler_ModuleI(t,x,p),tspan,x01,options);
t0 = [t1 ;t2+ t1(end)];
s0 = [s1;s2];
% 
% p(7) = p(7)*2;
% [t3,s3] = ode23s(@(t,x) ODE_Model_2PIControler_ModuleI(t,x,p),tspan,x01,options);
% t0 = [t0 ;t3+ t0(end)];
% s0 = [s0;s3];

figure(1)
subplot(3,1,1)
plot(t0([1,end]),[r2 r2],'-','Color',[200 200 200]*0.7/255,'LineWidth',3)
hold on
h2 = plot(t0,s0(:,1),'r',t0,s0(:,2),'b',t0,s0(:,3),'k','LineWidth',2)
legend([h2],'u_1','r_y','r_r'),hold off

subplot(3,1,2)
plot(t0([1,end]),[r1 r1],'-','Color',[200 200 200]/255,'LineWidth',3),hold on
h2 = plot(t0,s0(:,4),'r',t0,s0(:,5),'b',t0,s0(:,6),'k','LineWidth',2)
legend([h2],'u_2','r_','r_r2'),hold off

subplot(3,1,3)
h2 = plot(t0,s0(:,7),'r',t0,s0(:,8),'b',t0,s0(:,9),'k','LineWidth',2)
legend([h2],'y','w','z'),hold off



%%
hFig=figure(10);
set(hFig,'Units','inches', 'Position', [0 100 6 3]*1.3)

subplot(2,3,1:3)
plot(t0([1,end]),[r1 r1],'--','Color',[236 24 72]/255,'LineWidth',3)
hold on
plot(t0,s0(:,4),'Color',[144 189 49]/255,'LineWidth',3)
text(14, 0.27,strcat('2',var(7)),'FontSize',12);
text(29, 0.27,strcat('4',var(7)),'FontSize',12);
hold off
xlabel('time (h)','interpreter','latex')
ylabel('$y$ $(\mu M)$','interpreter','latex')
ax = gca;
ax.XTick = [0 tf  2*tf  3*tf ];
ax.YTick = [0 r1 yt];
xlim([0 tf*3]),ylim([0 yt])
N = 100000;
% e_v = linspace(0,1,N)*et;
u_v = logspace(-6,0,N)*ut; %e_v(end) = 0.99999*et;
y_v = logspace(-6,0,N)*yt;

y_e = A_F_u_I(u_v,[p1 0 r1]);
u_e = u_F_y(y_v,p2);
subplot(2,3,4),
plot(y_e,u_v,'Color',[153 153 153]/255,'LineWidth',4), hold on
plot(y_v,u_e,'Color',[244 128 36]/255,'LineWidth',4)
plot(s1(:,4),s1(:,1),'Color',[144 189 49]/255,'LineWidth',3)
hold off
xlabel('$y$ $(\mu M)$','interpreter','latex')
ylabel('$u$ $(\mu M)$','interpreter','latex')
axis([0 yt 0 ut])
ax = gca;
ax.XTick = [0  r1 yt];
ax.YTick = [0 round(s1(end,1)*10)/10 ut];

y_e = A_F_u_I(u_v,[p1 0 r1]);
p2(1)=2*p2(1);
u_e = u_F_y(y_v,p2);
subplot(2,3,5),
plot(y_e,u_v,'Color',[153 153 153]/255,'LineWidth',4), hold on
plot(y_v,u_e,'Color',[244 128 36]/255,'LineWidth',4)
plot(s2(:,4),s2(:,1),'Color',[144 189 49]/255,'LineWidth',3)
hold off
xlabel('$y$ $(\mu M)$','interpreter','latex')
ylabel('$u$ $(\mu M)$','interpreter','latex')
axis([0 yt 0 ut])
ax = gca;
ax.XTick = [0  r1 yt];
ax.YTick = [0 round(s2(end,1)*10)/10 ut];

y_e = A_F_u_I(u_v,[p1 0 r1]);
p2(1)=2*p2(1);
u_e = u_F_y(y_v,p2);
subplot(2,3,6),
plot(y_e,u_v,'Color',[153 153 153]/255,'LineWidth',4), hold on
plot(y_v,u_e,'Color',[244 128 36]/255,'LineWidth',4)
plot(s3(:,4),s3(:,1),'Color',[144 189 49]/255,'LineWidth',3)
hold off
xlabel('$y$ $(\mu M)$','interpreter','latex')
ylabel('$u$ $(\mu M)$','interpreter','latex')
axis([0 yt 0 ut])
ax = gca;
ax.XTick = [0  r1 yt];
ax.YTick = [0 round(s3(end,1)*10)/10 ut];
%%
% p1 = [k b a ph rh ut]; % 6 Controller 
% p2 = [k th g psi ph g1t A1t]; % 7 Plan
if(1)
%     gr1 =0.4;
    p = [p1 p2 r1];
    tf= 15;
tspan =[0 tf];


% Parameter Scan
vec1 = [1:6]; % 9
vec2 = [7:14]; % 9
Scale = log(2);
N = 9;
DIM1 = [0 6.5 10 4];
DIM2 = [0 0 12 4];
hFig=figure(3);
set(hFig,'Units','inches', 'Position', DIM1)
j = 1;
for i=vec1
    subplot(2,3,j)
    plot([0 tf],[r1 r1],'-','Color',[236 24 72]/255,'LineWidth',3)
    hold on,Plot_ScanParameter(p,i,Scale,N,tf,x0)
    title(var(i))
    j = j+1;
end
hold off

hFig=figure(4);
set(hFig,'Units','inches', 'Position', DIM2)
j = 1;
for i=vec2
    subplot(2,4,j)
    plot([0 tf],[r1 r1],'-','Color',[236 24 72]/255,'LineWidth',3)
    hold on ,Plot_ScanParameter(p,i,Scale,N,tf,x0)
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
    subplot(2,xx,j),plot([0 tf],[r1 r1],'-','Color',[153 153 153]/255,'LineWidth',3)
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
    subplot(2,xx,j),plot([0 tf],[r1 r1],'-','Color',[153 153 153]/255,'LineWidth',3)
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
    subplot(2,xx,j),plot([0 tf],[r1 r1],'-','Color',[236 24 72]/255,'LineWidth',3)
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
    subplot(2,xx,j),plot([0 tf],[r1 r1],'-','Color',[236 24 72]/255,'LineWidth',3)
    title(var(i))
    hold on ,Plot_ScanParameter(p,i,Scale,N,tf,x0,COLOR)
    subplot(2,xx,j+xx),Plot_ScanParameterNCC(p,i,Scale,N,COLOR)
    title(var(i))
    j = j+1;
end
hold off
