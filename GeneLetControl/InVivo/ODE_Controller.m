function y = ODE_Controller(x,p,A,I)
k = p(1);
b = p(2);
a = p(3);
ph = p(4);
g = p(5);
ut = p(6);

u = x(1);RA = x(2);RI = x(3);

cu = ut - u;

y(1) = a*RA*cu - b*RI*u - ph*u;
y(2) = k*A - a*RA*cu - g*RA*RI - ph*RA ;
y(3) = k*I - b*RI*u - g*RA*RI- ph*RI;
y = y(:);
end