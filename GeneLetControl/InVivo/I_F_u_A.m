function I_v = I_F_u_A(u,p)

% Kinetic parameters
k1 = p(1); % K1 = Kcat*Et
b = p(2);
a = p(3);
ph = p(4);
g = p(5);
ut = p(6);
A = p(7);
I = p(8);

A1 = b*g*u;
B1 = b*u.*(a*(ut-u)+ph) + ph*g*u;
C1 = -a*(ut-u)*k1*A + ph*u*(a*(ut-u)+ph);

rI = (-B1 +(B1.^2-4*A1.*C1).^0.5)./(2*A1);

rA = (b*u.*rI + ph*u)/a./(ut-u);

I_v = A - ph/k1*(rA-rI) - ph*u/k1;