function out = A_F_u_I(u1,u2,p)
k = p(1);
th = p(2);
% g = p(3);
% d = p(4);
kAp = p(5);
kAn = p(6);
kRp = p(7);
kRn = p(8);
n = p(9);
ut = p(10);
KA = p(11);
KI = p(12);

xA = Compute_Sequestration(p,k*u1./(u1+KA),th*u2./(u2+KI));
xR = Compute_Sequestration(p,th*u2./(u2+KI),k*u1./(u1+KA));

kA = kAn/kAp;
kR = kRn/kRp;

out = (xA.^n./(xA.^n/kA + xR.^n/kR + 1))/kA*ut;

end

function out = Compute_Sequestration(p,u1,u2)
g = p(3);
d = p(4);

B = (u2-u1)/d + d/g;
C = -u1/g;

out = (-B + (B.^2 - 4*C).^0.5)/2;
end

