% Solve RNA clock
%close all,clc, clear all
% Bistable parametrs
h = 3600;
uM = 10^(-6);
uMh = uM*h;
F = 1;
% Kinetic rate of biocontroller
kc = 20.*10^(-4)*h/4*F;
thc = 20.*10^(-4)*h/4*1.5*F;
bc = 3*10^(4)*uMh*2*1;
ac = 3*10^(4)*uMh*2*1;
phc = log(2)*60/30; 
gc = 3*10^4*uMh;
ut = 0.5*1;

pI = [kc thc bc ac phc gc ut]; % 5-11

% Kinetic rate of biocontroller plant
as = 20.*10^(-4)*h*.5/1;
ks = .2/1; % /2
phs = log(2)*60/3/1;
trs = .23*60/4*1; %*2
ds = log(2)*60/30;
ns = 2;
yt = 0.6;
pS = [as ks phs trs ds ns]; % 12-17

r1 = 0.1/1.5*2; % 2 3 4
R1 = thc*r1/kc;



c1 = 1;% proportional 
c2 = 1;% Integral
p = [pI pS r1];

Omega = 600/1;

set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultAxesFontSize', 12)

op = 1;
% Volume


% Get reactions
[ S, h, endSim ] = Model_Feedback( Omega,p );

% Initial condition
x0 = round(Omega*[0. 0.  0 0 0 ]')*0;
% x0 = (Omega*[0.1 0 0 0 0 0  0 0]');
% Simulation time
Tmax = 20;
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
% Write the initial state to the first "record"
X(:,1) = x;
T(1) = t;

[ S, h, endSim ] = Model_Feedback( Omega,p );
if (op==1)
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
else
step = 1;
while t<Tmax
    rates = h(x);
    x = x+S*rates*deltaT;
    t = t+deltaT;
    step = step + 1;
    X(:,step) = x;
    T(step) = t;
end
% Deterministic
end

 X = X/Omega;
% Plot the trajectory accumulated in arrays T and X:
% r1(1) - r4(2) - g2(3) - cr1(4) - r2(5) - r3(6) - g1(7) - cr3(8)
%
set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultAxesFontSize', 12)

figure(2)
subplot(3,1,1)
plot(T,X(2,:),'Color',[255 59 48]/255,'LineWidth',3),hold on
plot(T,X(3,:),'Color',[0 122 255]/255,'LineWidth',3),hold off
legend('r_y','r_r')
xlabel('time (h)')
ylabel('\mu M')
xlim([0 Tmax])

subplot(3,1,2)
% plot(T,X(3,:),'Color',[255 59 48]/255,'LineWidth',3),hold on
plot(T,X(1,:),'Color',[0 122 255]/255,'LineWidth',3),hold off
legend('u')
xlabel('time (h)')
ylabel('\mu M')
xlim([0 Tmax])

subplot(3,1,3)
plot([0 Tmax],[R1 R1],'k-','LineWidth',1),hold on
plot(T,X(5,:),'Color',[0 122 255]/255,'LineWidth',3),hold off
legend('r','y')
xlabel('time (h)')
ylabel('\mu M')
xlim([0 Tmax])

%%
figure(22)
subplot(1,2,1)
plot([0 Tmax],[R1 R1],'k-','LineWidth',1),hold on
plot(T,X(5,:),'Color',[0 122 255]/255,'LineWidth',2),hold off
legend('r','y')
xlabel('Time (h)')
ylabel('\mu M')
xlim([0 Tmax])
ylim([0 0.6])
