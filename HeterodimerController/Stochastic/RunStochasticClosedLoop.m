% clear all,clc
%///////////////////////////////
% Volume
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

r1 = 0.2; % 0.1/2/4
R1 = KIc/(thc/kc*(1+KAc/r1)-1); 

p1 = [kc thc gc phc bc kAp kAn kRp kRn nc ut KAc KIc]; % 1 - 11


% Kinetic rate of biocontroller plant
as = 20.*10^(-4)*h*.5;
phs = log(2)*60/3;
trs = .23*60/4;
ds = log(2)*60/30;


p2 = [as phs trs ds]; % 12-15

set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultAxesFontSize', 12)

p = [p1 p2 r1];
Omega = 600;
% Get reactions
[ S, h, endSim ] = Model_ClosedLoop( Omega,p );

% Initial condition
x0 = round(Omega*[0. 0 0 0.  0. 0]');

% Simulation time
Tmax = 30;
%///////////////////////////////


% Suggested time step
deltaT = 0.01;

% Initialize simulation:
x = x0;
t = 0;


% Initialize "recording" variables to save the state:
% (It's much faster to reserve the whole memory at this point
% instead of enlarging the array in every step)
X = nan(size(x,1),100000);
T = nan(1,100000);
u =nan(1,100000);
% Write the initial state to the first "record"
X(:,1) = x;
T(1) = t;
u(1) = 0;

step = 1;
while t<Tmax
    hs = h(x);
    hscs = cumsum(hs);
    rateall = sum(hs);
    tau = -log(rand)/rateall;
    r = rand*sum(hs);
    j = find(r<=hscs,1,'first');
    x = x + S(:,j);
    t = t+tau;
    step = step + 1;
    X(:,step) = x;
    T(step) = t;
end


% Plot the trajectory accumulated in arrays T and X:
% r1(1) - r4(2) - g2(3) - cr1(4) - r2(5) - r3(6) - g1(7) - cr3(8)
%
% figure(1);
% 
% h1a = plot(T,X(1,:),'Color',[255 110 30]/255,'LineWidth',3), hold on
% h2a = plot(T,X(2,:),'Color',[255 110 30]/255,'LineWidth',3);h2a.Color(4)=0.4;
% plot(T,X(7,:),'Color',[1 1 1]*0.5,'LineWidth',3);h2a.Color(4)=0.4;
% hold off
% legend('x_1','x_2','u')
% % legend('x_1','u')
% xlabel('time (h)'), ylabel('[\mu M]')
% xlim([0 Tmax])
 %% Plot 2
 
 X = X/Omega;



%%
%
NN = 100000;
y_c = linspace(0,1,NN)*1; %e_v(end) = 0.99999*et;
u_c = A_F_u_I(r1,y_c,p1);
idx = find (y_c<0);y_c(idx) = [];y_c(idx)=[];

u_s = logspace(-6,0,NN);
y_s = u_F_y(u_s,p2);

hFig=figure(3);
set(hFig,'Units','inches', 'Position', [0 100 7/2 2*1])
subplot(1,2,1)
plot([0 Tmax],[R1 R1],'k-','LineWidth',1), hold on
plot(T,X(6,:),'-','Color',[113 111 178]/255,'LineWidth',2)
xlabel('time(h)'),ylabel('y (\mu M)')
hold off
subplot(1,2,2)
plot(y_s,u_s,'Color',[160 158 158]/255,'LineWidth',2), hold on
plot(y_c,u_c,'Color',[248 152 56]/255,'LineWidth',2)
plot(X(6,:),X(3,:),'Color',[113 111 178]/255,'LineWidth',2)
axis([0 1 0 ut])
hold off
xlabel('y (\mu M)'),ylabel('u (\mu M)')
% axis([0 0.05 0 .5])