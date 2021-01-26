function u_v = u_F_y(y,p)

% Kinetic parameters
a = p(1);
ku = p(2);
ph = p(3);
tr = p(4);
d = p(5);
n = p(6);

ry = d*y/tr;
u_v = ku.*(a./(ph*ry) - 1).^(1/n);
% u_v = d*pu/tr;
% u_v = (b*y + ph + g*w).*z/k;