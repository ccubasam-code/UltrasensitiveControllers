function dx = ODE_Gene(x,p,u)
% The input will inhibited the gene "y"
a = p(1);
ku = p(2);
ph = p(3);
tr = p(4);
d = p(5);
n = p(6);


ry = x(1); y = x(2);

dx(1) = a*ku^n/(u^n+ku^n) - ph*ry;
dx(2) = tr*ry - d*y;
dx = dx(:);
end