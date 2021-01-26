function [w, z] = Compute_EQ_Proc(y,p)

% Kinetic parameters
k = p(1);
b = p(2);
a = p(3);
ph = p(4);
g = p(5);
th = p(6);
yt = p(7);
wt = p(8);

cy = yt-y;

A1 = a*g*cy/b./y;
B1 = a*cy + th;
C1 = -th*(wt-y);

w = (-B1 +(B1.^2-4*A1.*C1).^0.5)./(2*A1);

z = a*cy/b./y.*w;
