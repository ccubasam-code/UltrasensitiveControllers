function out = A_F_u_I(u1,u2,p)
k = p(1);
th = p(2);
d = p(4);
b = p(5);
kAp = p(6);
kAn = p(7);
kRp = p(8);
kRn = p(9);
n = p(10);
ut = p(11);
KA = p(12);
KI = p(13);

kA = kAn/kAp;
kR = kRn/kRp;

xA = Compute_Sequestration(p,k*u1./(u1+KA),th*u2./(u2+KI));
xR = Compute_Sequestration(p,th*u2./(u2+KI),k*u1./(u1+KA));
out = (xA.^n./(xA.^n/kA + xR.^n/kR + 1))/kA*ut;

end

function out = Compute_Sequestration(p,u1,u2)
g = p(3);
d = p(4);

B = (u2-u1)/d + d/g;
C = -u1/g;

out = (-B + (B.^2 - 4*C).^0.5)/2;
end

