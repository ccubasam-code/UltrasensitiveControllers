h = 3600;
uM = 10^(-6);
uMh = uM*h;
% Kinetic rate of biocontroller
kc = 20.*10^(-4)*h/4;
thc = 20.*10^(-4)*h/4*1.5;
bc = 3*10^(4)*uMh*2*2;
ac = 3*10^(4)*uMh*2*2;
phc = log(2)*60/30; 
gc = 3*10^4*uMh;
ut = 0.5*1;

trc = .23*60/4;
dc = log(2)*60/30; % 30 minutes

r1 = 0.1/1.5;
r2 = 0.2/1.5;
r3 = 0.3/1.5;
R1 = thc*r1/kc;
R2 = thc*r2/kc;
R3 = thc*r3/kc;

pI = [kc thc bc ac phc gc ut 0 r1]; % 1-7
pC = [kc thc phc gc 0 r1];

as = 20.*10^(-4)*h*.5/1*2;
ks = .2/1; % /2
phs = log(2)*60/3;
trs = .23*60/4*1; %*2
ds = log(2)*60/30;
ns = 2;

pS = [as ks phs trs ds ns]; % 8-15

m = 0.5*bc/phc*ut/1;

N = 1000;
y_S = logspace(-3,0,N)*trs/ds*as/phs;
x_S = u_F_y(y_S,pS);

N = 100000;
% e_v = linspace(0,1,N)*et;
u_v = logspace(-6,0,N)*ut; %e_v(end) = 0.99999*et;

% Tracking
%
x = linspace(0,1,N);
y_I = ut*x.^m./(R1.^m+x.^m);
y_P = Func_IO_Titration(x,pC);

range = [0 0.6 0 0.8];

hFig=figure(11);
set(hFig,'Units','inches', 'Position', [0 100 3.5*2 2])
Col1 = [160 158 158]/255;
Col2 = [248 152 56]/255;
subplot(1,3,1)
plot(x_S,y_S,'Color',Col1,'LineWidth',2),hold on
plot(x,y_P,'Color',Col2,'LineWidth',2),hold off
axis(range)
ax = gca;ax.XTick = [0 R1  0.6 ];ax.YTick = [0 0.8 ];
xlabel('Time')
subplot(1,3,2)
plot(x_S,y_S,'Color',Col1,'LineWidth',2),hold on
plot(x,y_I,'Color',Col2,'LineWidth',2),hold off
axis(range)
ax = gca;ax.XTick = [0 R1  0.6 ];ax.YTick = [0 0.8 ];
xlabel('Time')
subplot(1,3,3)
plot(x_S,y_S,'Color',Col1,'LineWidth',2),hold on
plot(x,y_P+y_I,'Color',Col2,'LineWidth',2),hold off
axis(range)
ax = gca;ax.XTick = [0 R1  0.6 ];ax.YTick = [0 0.8 ];
xlabel('Time')