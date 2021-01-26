function I_v = I_F_u_A(u,p)

% Kinetic parameters
k = p(1);
th = p(2);
b = p(3);
a = p(4);
ph = p(5);
g = p(6);
ut = p(7);
A = p(8);
I = p(9);

A1 = b*g*u;
B1 = b*u.*(a*(ut-u)+ph);
C1 = -a*(ut-u)*k*A;

rI = (-B1 +(B1.^2-4*A1.*C1).^0.5)./(2*A1);

rA = b*u.*rI/a./(ut-u);

I_v = k*A/th - ph/th*(rA-rI);