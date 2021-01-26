function y_v = u_F_y(u,p)

% Kinetic parameters
a = p(1);
ph = p(2);
psi = p(3);
d = p(4);


y_v = a/ph*psi/d*u;
