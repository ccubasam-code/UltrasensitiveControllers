function dx = ODE_GeneletI(x,p,u)
% The input will inhibited the gene "y"
k = p(1);
b = p(2);
a = p(3);
ph = p(4);
g = p(5);
th = p(6);
yt = p(7);
wt = p(8);

y = x(1); w = x(2); z = x(3);

cy = yt - y;
cw = wt - w - y;

dx(1) = a*w*cy - b*y*z;
dx(2) = th*cw - a*cy*w - g*w*z;
dx(3) = k*u - b*y*z - ph*z - g*w*z;
dx = dx(:);
end