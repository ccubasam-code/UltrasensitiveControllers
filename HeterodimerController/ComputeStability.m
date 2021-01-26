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

function [J,tmp] = ComputeJacobian(uA,y,p)
th = p(2);
gc = p(3);
dc = p(4);
kAp = p(5);
kAn = p(6);
kRp = p(7);
kRn = p(8);
ut = p(10);
KI = p(12);

as = p(13);
phs = p(14);
psi = p(15);
ds = p(16);

r = p(17);

[xA, xR] = ComputeEqConntroller(r,y,p);
[m] = ComputeEqProcess(y,p);

kA = kAn/kAp;
kR = kRn/kRp;

uR = (xR^2*kA/(xA^2*kR))*uA;

u = ut - uA - uR; 

tmp = [xA xR uA uR m y];


J = zeros(6,6);

J(1,1) = -(gc*xR + dc) - 4*kAp*xA*u;
J(1,2) = -gc*xA;
J(1,3) = 2*(kAp*xA^2 + kAn);
J(1,4) = 2*kAp*xA^2;

J(2,1) = -gc*xR;
J(2,2) = -(gc*xA + dc) - 4*kRp*xR*u;
J(2,3) = 2*kRp*xR^2;
J(2,4) = 2*(kRp*xR^2 + kRn);
J(2,6) = th*KI/(y + KI)^2;

J(3,1) = 2*kAp*xA*u;
J(3,3) = -kAp*xA^2 - kAn;
J(3,4) = -kAp*xA^2;

J(4,2) = 2*kRp*xR*u;
J(4,3) = -kRp*xR^2;
J(4,4) = -kRp*xR^2 - kRn;

J(5,3) = as;
J(5,5) = -phs;

J(6,5) = psi;
J(6,6) = -ds;

end

function [xA, xR] = ComputeEqConntroller(u1,u2,p)
k = p(1);
th = p(2);

KA = p(11);
KI = p(12);


xA = Compute_Sequestration(p,k*u1./(u1+KA),th*u2./(u2+KI));
xR = Compute_Sequestration(p,th*u2./(u2+KI),k*u1./(u1+KA));

end

function out = Compute_Sequestration(p,u1,u2)
g = p(3);
d = p(4);

B = (u2-u1)/d + d/g;
C = -u1/g;

out = (-B + (B.^2 - 4*C).^0.5)/2;
end

function [m] = ComputeEqProcess(y,p)
psi = p(15);
d = p(16);

m = d*y/psi;
end