function A_v = A_F_u_I_Break(u,p)
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
w = p(12);

ut = rh/d;

d = d*op;


rA = (b*w*u+ d*u)/a./(ut-u);
rI = th*I./(g*rA+ph);

A_v = th*I/k + ph*(rA-rI)/k + a*rA.*(ut-u)/k;
% A_v = th*I/k + ph*(rA-rI)/k + (b*w*u+d*u)/k;