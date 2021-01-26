function [rA, rI] = Compute_EQ_Contr(u,p)

% Kinetic parameters
k = p(1);
th = p(2);
b = p(3);
a = p(4);
ph = p(5);
g = p(6);
ut = p(7);
% A = p(8);
I = p(9);

A1 = a*g*(ut-u);
B1 = a*(ut-u).*(b*u+ph);
C1 = -b*u*th*I;

rA = (-B1 +(B1.^2-4*A1.*C1).^0.5)./(2*A1);
rI = a*(ut-u).*rA./(b*u);
% A_v = th*I/k + ph*(rA-rI)/k;