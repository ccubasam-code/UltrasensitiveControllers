function dx = ODE_PI_ClosedLoop(t,x,p)
% Return T - RI - RA
% Parameters
pP = p(1:4);
pI = p(5:12);
pS = p(13:18);

c1 = p(19); % Reference
c2 = p(20);
r = p(21);

dx = nan(1,7);
% States
y = x(7);
% u - RA - RI
dx(1:2) = ODE_Proportional(x(1:2),pP,y,r);
dx(3:5) = ODE_ControllerDycU(x(3:5),pI,y,r);
u1 = x(1);
u2 = x(3);
% y - w - z
dx(6:7) = ODE_Gene(x(6:7),pS,u1*c1+u2*c2);
dx = dx(:);
end
