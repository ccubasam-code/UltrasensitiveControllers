% Solve RNA clock
%close all,clc, clear all
% Bistable parametrs
h = 3600;
uM = 10^(-6);
uMh = uM*h;
% Kinetic rate of biocontroller
kc = 20.*10^(-4)*h/4/10;
thc = 20.*10^(-4)*h/4/10;
bc = 3*10^6*uMh*2;
ac = 3*10^6*uMh*2;
phc = log(2)*60/3; % 30 minutes
gc = 3*10^4*uMh;
ut = 0.5;
KAc = 0.4;
KIc = 0.3;
trc = .23*60/4;
dc = log(2)*60/30; % 30 minutes

r1 = 0.15;
y0 = 0.017;
% y0 = 0.5;
R1 = KAc/(kc/thc*(1+KIc/r1)-1); 
yt=1;

% p1 = [kc thc bc ac phc gc ut KAc KIc trc dc]; % 1-11
p1 = [kc thc bc ac phc gc ut KAc KIc]; % 1-9


% Kinetic rate of biocontroller plant
as = 20.*10^(-4)*h*.5;
ks = .2;
phs = log(2)*60/3;
trs = .23*60/4;
ds = log(2)*60/30;
ns = 2*1;

p2 = [as ks phs trs ds ns]; % 12-17

var = {'\kappa_c','\beta_c','\alpha_c','\phi_c','\gamma_c', 'u^{tot}', ...
    '\kappa_s','\beta_s','\alpha_s','\phi_s','\gamma_s','\theta_s','y^{tot}','w^{tot}'};
set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultAxesFontSize', 12)

p = [p1 p2 r1];


%
N = 100000;
% e_v = linspace(0,1,N)*et;
u_v = logspace(-6,0,N)*ut; %e_v(end) = 0.99999*et;
y_v = logspace(-6,0,N)*yt;

y_e = A_F_u_I(u_v,[p1 0 r1]); % controller
u_e = u_F_y(y_v,p2); % Process

idx=find(y_e>yt);
y_e = y_e(1:idx(1));
u_v = u_v(1:idx(1));
% Iteraction map
nn = 10;

M = nan(nn,2); % y - u
M(1,:) = [y0 0];
for i=2:2:nn
    M(i,:) = [M(i-1,1) u_F_y(M(i-1,1),p2)]; % Process
    M(i+1,:) = [A_F_u_I(M(i,2),[p1 0 r1]) M(i,2)]; % Process
end


hFig=figure(2);
set(gcf,'PaperPositionMode','auto')
set(hFig,'Units','inches', 'Position', [0 10 7/2 2])

plot(y_e,u_v,'Color',[248 152 56]/255,'LineWidth',2), hold on
plot(y_v,u_e,'Color',[160 158 158]/255,'LineWidth',2)
plot(M(:,1),M(:,2),'b-','LineWidth',1)
hold off
xlabel('$y$ $(\mu M)$','interpreter','latex')
ylabel('$u$ $(\mu M)$','interpreter','latex')
% title (['|G_cG_P| = ', num2str(round(10*G)/10)])
axis([0 0.6 0 ut])
ax = gca;
ax.XTick = [0  R1 0.6];
% ax.YTick = [0 round(s1(end,1)*10)/10 ut];
ax.YTick = [0  ut];

%

hFig=figure(3);
set(gcf,'PaperPositionMode','auto')
set(hFig,'Units','inches', 'Position', [0 10 7/2 2])

plot(y_e,u_v,'Color',[248 152 56]/255,'LineWidth',2), hold on
plot(y_v,u_e,'Color',[160 158 158]/255,'LineWidth',2)
xlabel('$y$ $(\mu M)$','interpreter','latex')
ylabel('$u$ $(\mu M)$','interpreter','latex')
% title (['|G_cG_P| = ', num2str(round(10*G)/10)])
axis([0 0.6 0 ut])
ax = gca;
ax.XTick = [0  R1 0.6];
% ax.YTick = [0 round(s1(end,1)*10)/10 ut];
ax.YTick = [0  ut];
plot(M(:,1),M(:,2),'b-','LineWidth',1)

y0 = 0.1;
M = nan(nn,2); % y - u
M(1,:) = [y0 0];
for i=2:2:nn
    M(i,:) = [M(i-1,1) u_F_y(M(i-1,1),p2)]; % Process
    M(i+1,:) = [A_F_u_I(M(i,2),[p1 0 r1]) M(i,2)]; % Process
end
plot(M(:,1),M(:,2),'r-','LineWidth',1)

y0 = 0.5;
M = nan(nn,2); % y - u
M(1,:) = [y0 0];
for i=2:2:nn
    M(i,:) = [M(i-1,1) u_F_y(M(i-1,1),p2)]; % Process
    M(i+1,:) = [A_F_u_I(M(i,2),[p1 0 r1]) M(i,2)]; % Process
end
plot(M(:,1),M(:,2),'k-','LineWidth',1)
hold off