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
thc = 20.*10^(-4)*h/4*1.1*I; % 1.5
bc = 3*10^4*uMh*2;
ac = 3*10^4*uMh*2;
phc = log(2)*60/30; % 3 minutes 
gc = 3*10^4*uMh;

ut = 0.5;




Omega = 600;

p = [kc thc bc ac phc gc ut]; % 6

% Get reactions
[ S, h, endSim ] = Model_Brink( Omega,p );

% Initial condition
x0 = round(Omega*[ 0 0 0 ]');

% Simulation time
Tmax = 5*1;
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


%% Plot 3
figure(50);
h1a = plot(T,X(1,:),'Color',[255 110 30]/255,'LineWidth',3), hold on
plot(T,X(2,:),'Color',[1 1 1]*0.5,'LineWidth',3);
hold off
legend('x_1','y_1')
xlabel('time (h)'), ylabel('\mu M')
xlim([0 Tmax])

figure(51);
h1a = plot(T,X(1,:),'Color',[255 110 30]/255,'LineWidth',3), hold on
plot(T,X(2,:),'Color',[1 1 1]*0.5,'LineWidth',3);
plot(T,X(3,:),'-.','Color',[1 1 1]*0.3,'LineWidth',3);
% plot(T,X(4,:),'r','LineWidth',3);
hold off
legend('x_1','y_1','x_2','y_2')
xlabel('time (h)'), ylabel('\mu M')
xlim([0 Tmax])

% ylim([0 2])
figure(6);
h1a = plot(T,X(1,:),'Color',[255 110 30]/255,'LineWidth',3), hold on
h2a = plot(T,X(3,:),'Color',[255 110 30]/255,'LineWidth',3);h2a.Color(4)=0.4;
% plot(T,X(5,:),'Color',[1 1 1]*0.5,'LineWidth',3);
hold off
legend('x_1','x_2')
xlabel('time (h)'), ylabel('\mu M')
xlim([0 Tmax])
%
figure(7);
h1a = plot(T,X(2,:),'Color',[255 110 30]/255,'LineWidth',3), hold on
h2a = plot(T,X(3,:),'Color',[1 1 1]*0.5,'LineWidth',3);h2a.Color(4)=1;
hold off
legend('y_1','y_2')
xlabel('time (h)'), ylabel('\mu M')
xlim([0 Tmax])

figure(8);
h1a = plot(X(1,:),X(3,:),'Color',[255 110 30]/255,'LineWidth',3)
% legend('x_1+x_2','u')
xlabel('x_1'), ylabel('x_2')
% xlim([0 Tmax])