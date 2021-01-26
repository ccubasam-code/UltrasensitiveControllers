h = 3600;
uM = 10^(-6);
uMh = uM*h;

P1 = 0.01;
P2 = 0.01;
lam = [0 1 2];
g = 5*10^9*uMh;
ph = log(2)/30*60;

p = [P1 P2 lam(1) g ph];

u1 = .1;
u2 = linspace(0,1,501);u2(1)=[];

xn0 = Compute_IOV2(p,u1,u2);
figure(1)
subplot(2,1,1)
plot(u2,xn0,'b-','LineWidth',2)

hold on
p = [P1 P2 lam(2) g ph];
xn1 = Compute_IOV2(p,u1,u2);
plot(u2,xn1,'r','LineWidth',2)

p = [P1 P2 lam(3) g ph];
xn2 = Compute_IOV2(p,u1,u2);
plot(u2,xn2,'k','LineWidth',2)
hold off
xlim([0  2*u1])


subplot(2,1,2),plot(u2/u1,xn0/max(xn0),'b-','LineWidth',2)
hold on
plot(u2/u1,xn1/max(xn1),'r-','LineWidth',2)
plot(u2/u1,xn2/max(xn2),'k-','LineWidth',2)
hold off
xlim([0  2])


%
u2 = .1;
u1 = linspace(0,1,501);u1(1)=[];
p = [P1 P2 lam(1) g ph];
xn0 = Compute_IOV2(p,u1,u2);
figure(2)
subplot(2,1,1)
plot(u1,xn0,'b-','LineWidth',2)

hold on
p = [P1 P2 lam(2) g ph];
xn1 = Compute_IOV2(p,u1,u2);
plot(u1,xn1,'r','LineWidth',2)

p = [P1 P2 lam(3) g ph];
xn2 = Compute_IOV2(p,u1,u2);
plot(u1,xn2,'k','LineWidth',2)
hold off
xlim([0  2*u2])

subplot(2,1,2)
plot(u1/u2,xn0/max(xn0),'b-','LineWidth',2)
hold on
plot(u1/u2,xn1/max(xn1),'r-','LineWidth',2)
plot(u1/u2,xn2/max(xn2),'k-','LineWidth',2)
hold off
xlim([0  2])