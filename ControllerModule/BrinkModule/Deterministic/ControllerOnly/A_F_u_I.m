function A_v = A_F_u_I(u,p)
% This consider production and degradation of species u

% Kinetic parameters
k = p(1);
th = p(2);
b = p(3);
a = p(4);
ph = p(5);
g = p(6);
rh = p(7);
d = p(8);
A = p(9);
I = p(10);
op = p(11);


ut = rh/d;
% ut = 0.5;

d = d*op;

A1 = a*g*(ut-u);
B1 = a*(ut-u).*(b*u+ph) - d*u*g;
C1 = -b*u*th*I - d*u.*(b*u+ph);

rA = (-B1 +(B1.^2-4*A1.*C1).^0.5)./(2*A1);

rI = (a*(ut-u).*rA - d*u)./(b*u);

A_v = th*I/k + ph*(rA-rI)/k + d/k*u;