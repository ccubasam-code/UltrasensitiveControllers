function y = ODE_ControllerSeq(x,p,A,I)
k = p(1);
th = p(2);
b = p(3);
a = p(4);
ph = p(5);
g = p(6);
ut = p(7);
z = 0.1/10; % replace the expression of RI 
u = x(1);RA = x(2);RI = x(3);

cu = ut - u;

y(1) = a*RA*cu - b*z*u;
y(2) = k*A - a*RA*cu - g*RA*RI - ph*RA*0 ;
y(3) = th*I - b*RI*u*0 - g*RA*RI- ph*RI*0;
y = y(:);
end