% Solve RNA clock
%close all,clc, clear all
% Bistable parametrs
h = 3600;
uM = 10^(-6);
uMh = uM*h;
set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultAxesFontSize', 10)
% Kinetic rate of biocontroller
A = 0.1;
I = 0.1;
ii = 1;

kc = 20.*10^(-4)*h/4*A;
thc = 20.*10^(-4)*h/4*1.5*I;
bc = 3*10^4*uMh*2;
ac = 3*10^4*uMh*2;
phc = log(2)*60/30; % 3 minutes 
gc = 3*10^4*uMh;

ut = 0.5;




Omega = 600;

p = [kc thc bc ac phc gc ut]; % 6
p0 = p;
vec = [0.5 0.8 0.9 1 1.1 1.2 2]*1.5;
for ii =1:7
    p0(1) = p(1)*vec(ii);
% Get reactions
[ S, h, endSim ] = Model_Brink( Omega,p0 );

% Initial condition
x0 = round(Omega*[ 0 0 0 ]');

% Simulation time
Tmax = 5*4;
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


 X = X/Omega;
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

% Plot 2

hFig=figure(22);
set(hFig,'Units','inches', 'Position', [0 100 7 1.5])

subplot(1,7,ii)
h1a = plot(T,X(3,:),'Color',[255 110 30]/255,'LineWidth',1), hold on
hold off
% legend('x_1','y_1')
% xlabel('time (h)'), ylabel('\mu M')
xlim([0 Tmax])
ylim([0 0.61])

ax = gca;
ax.XTick = [0 Tmax ];
ax.YTick = [0 2 4 6]*0.1;

end
