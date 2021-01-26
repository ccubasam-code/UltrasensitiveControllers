function [ExistEquilibrium,eigJ,eq] = ComputeStability(p,u,y)

if ~isempty(u)
    [J,eq] = ComputeJacobian(u,y,p);
    eigJ = eig(J);
    ExistEquilibrium =1;
else
    eigJ = zeros(6,1);
    ExistEquilibrium = 0;
    eq = zeros(6,1);
end

end

function [J,tmp] = ComputeJacobian(u,y,p)

kc = p(1);
thc = p(2);
bc = p(3);
ac = p(4);
phc = p(5); 
gc = p(6);
ut = p(7);

ks = p(8);
bs = p(9);
as = p(10);
phs = p(11);
gs = p(12);
ths = p(13);
yt = p(14);
wt = p(15);

[rY, rR] = ComputeEqConntroller(u,p,p(16));
[w,z] = ComputeEqProcess(y,p);
tmp = [u rY rR y w z];

% check jacobian
J = zeros(6,6);
J(1,1) = -(ac*(ut-u) + gc*rR + phc);
J(1,2) = gc*rY;
J(1,3) = ac*rY;
J(1,4) = kc;

J(2,1) = gc*rR;
J(2,2) = -(bc*u + gc*rY + phc);
J(2,3) = bc*rR;

J(3,1) = ac*(ut-u);
J(3,2) = bc*u;
J(3,3) = -(ac*rY+bc*rR);

J(4,4) = -as*yt*w/y;
J(4,5) = as*(yt-y);
J(4,6) = bs*y;

J(5,4) = as*w-ths;
J(5,5) = -ths*wt/w;
J(5,6) = gs*w;

J(6,3) = -ks;
J(6,4) = bs*z;
J(6,5) = gs*z;
J(6,6) = -ks*u/z;

end

function [rA, rI] = ComputeEqConntroller(u,p,r)
k = p(1);
th = p(2);
b = p(3);
a = p(4);
ph = p(5);
g = p(6);
ut = p(7);

A1 = a*g*(ut-u);
B1 = a*(ut-u).*(b*u+ph);
C1 = -b*u*th*r;

rA = (-B1 +(B1.^2-4*A1.*C1).^0.5)./(2*A1);
rI = a*(ut-u).*rA./(b*u);
end

function [w,z] = ComputeEqProcess(y,p)
ks = p(8);
bs = p(9);
as = p(10);
phs = p(11);
gs = p(12);
ths = p(13);
yt = p(14);
wt = p(15);

cy = yt-y;

A1 = as*gs*cy/bs./y;
B1 = as*cy + ths;
C1 = -ths*(wt-y);
w = (-B1 +(B1.^2-4*A1.*C1).^0.5)./(2*A1);
z = as*cy/bs./y.*w;
end