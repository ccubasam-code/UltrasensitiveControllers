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

r1 = .3/1.5;
R1 = thc*r1/kc;
pI = [kc thc bc ac phc gc ut 0 r1]; % 1-7
pC = [kc thc phc gc 0 r1];

as = 20.*10^(-4)*h*.5;
ks = .2;
phs = log(2)*60/3;
trs = .23*60/4;
ds = log(2)*60/30;
ns = 2;

p2 = [as ks phs trs ds ns]; % 8-15



N = 1000;
y_s = logspace(-3,0,N)*trs/ds*as/phs;
u_s = u_F_y(y_s,p2);


N = 100000;
% e_v = linspace(0,1,N)*et;
u_v = logspace(-6,0,N)*ut; %e_v(end) = 0.99999*et;

A_v = A_F_u_I(u_v,pI);

figure(1)
subplot(3,1,1)
plot(A_v,u_v,'k'), hold on

%
x = linspace(0,1,N);
m = 0.5*bc/phc*ut/1;
y_t = ut*x.^m./(R1.^m+x.^m);

plot(x,y_t,'b'), hold off
xlim([0 0.5])
ylim([0 ut])

%
subplot(3,1,2)

y_v = Func_IO_Titration(x,pC);
plot(x,y_v,'k')
xlim([0 0.5])
% ylim([0 ut])

subplot(3,1,3)

plot(x,y_v+y_t,'b'), hold on
plot(y_s,u_s,'k'), hold off
xlim([0 0.5])
ylim([0 ut*2])