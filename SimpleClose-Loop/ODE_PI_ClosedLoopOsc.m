function dx = ODE_PI_ClosedLoopOsc(t,x,p)
% Return T - RI - RA
% Parameters
pP = p(1:4);
pI = p(5:11);
pS = p(12:19);

c1 = p(20); % Reference
c2 = p(21);
r = p(22);

A = p(23);
T = p(24);
pS(1)=p(12) + A*(sin(2*pi/T*t+pi/2)+1)/2;
dx = nan(1,8);
% States
y = x(6);
% u - RA - RI
dx(1:2) = ODE_Proportional(x(1:2),pP,y,r);
dx(3:5) = ODE_Controller(x(3:5),pI,y,r);
u1 = x(1);
u2 = x(3);
% y - w - z
dx(6:8) = ODE_GeneletI(x(6:8),pS,u1*c1+u2*c2);
dx = dx(:);
end
