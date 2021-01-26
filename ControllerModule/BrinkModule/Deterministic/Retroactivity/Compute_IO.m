function xn = Compute_IO(p,u1,u2)

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
C = -a;

yd = (-B+(B.^2 -4*A.*C).^0.5)./(2*A);
zd = yd + (a-th)/ph;

f = P1*zd./(P2*yd);

L = 1 + lam;

AN = 1-f*L;
BN = -(f*(L*P2-1)-1-P1);
CN = f*P2;

xn = (-BN+(BN.^2 -4*AN.*CN).^0.5)./(2*AN);
% tmp = find(xn<0);
% xn(tmp)=tmp*0;

