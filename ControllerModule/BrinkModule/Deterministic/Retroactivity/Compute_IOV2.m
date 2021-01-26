function xn = Compute_IOV2(p,u1,u2)
% Improve computation of IO, V1 had negative values
P1 = p(1);
P2 = p(2);
lam = p(3);
g = p(4);
ph = p(5);

%%%%
a = u1;
th = u2;

A =  g;
B = g*(a-th)/ph + ph;
C = -th;

yd = (-B+(B.^2 -4*A.*C).^0.5)./(2*A);
zd = yd + (a-th)/ph;

% zd = a./(g*yd+ph);
f = P1*zd./(P2*yd);
% figure(10)
% subplot(2,1,1),plot(u2,yd,'r','LineWidth',2)
% subplot(2,1,2),plot(u2,zd,'r','LineWidth',2)
L = 1 + lam;



xn = f.^1./(f.^1*L+1);
% tmp = find(xn<0);
% xn(tmp)=tmp*0;

