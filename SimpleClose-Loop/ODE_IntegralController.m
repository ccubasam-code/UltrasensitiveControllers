function y = ODE_IntegralController(x,p,A,I)
k = p(1);
th = p(2);
ph = p(3);
g = p(4);

RA = x(1);RI = x(2);

%y(1) = a*RA*cu - b*RI*u;
y(1) = k*A  - g*RA*RI - ph*RA*0 ;
y(2) = th*I - g*RA*RI- ph*RI*0;
y = y(:);
end