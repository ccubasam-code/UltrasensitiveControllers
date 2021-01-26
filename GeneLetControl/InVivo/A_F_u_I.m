function A_v = A_F_u_I(u,p)

% Kinetic parameters
k1 = p(1); % K1 = Kcat*Et
b = p(2);
a = p(3);
ph = p(4);
g = p(5);
ut = p(6);
A = p(7);
I = p(8);

A1 = a*g*(ut-u);
B1 = a*(ut-u).*(b*u+ph) - g*ph*u;
C1 = -b*u*k1*I - ph*(b*u+ph).*u;

rA = (-B1 +(B1.^2-4*A1.*C1).^0.5)./(2*A1);

rI = (a*(ut-u).*rA - ph*u)./(b*u);
A_v = I + ph*(rA-rI)/k1 + ph*u/k1;