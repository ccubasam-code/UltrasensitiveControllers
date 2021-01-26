function y = ODE_ControllerDycU(x,p,A,I)
k = p(1);
th = p(2);
b = p(3);
a = p(4);
ph = p(5);
g = p(6);
rh = p(7);
d = p(8);

ut = rh/d;

u = x(1);RA = x(2);RI = x(3);

cu = ut - u;

y(1) = ut*ph*0 + a*RA*cu - b*RI*u - d*u;
y(2) = k*A - a*RA*cu - g*RA*RI - ph*RA ;
y(3) = th*I - b*RI*u - g*RA*RI- ph*RI;
y = y(:);
end