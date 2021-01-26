function dx = ODE_Model_PIControlerSeq_ModuleI(t,x,p)
% Return T - RI - RA
% Parameters
p1 = p(1:7);
p2 = p(8:15);
r = p(16); % Reference

dx = nan(1,6);
% States
y = x(4);
% u - RA - RI
dx(1:3) = ODE_ControllerSeq(x(1:3),p1,y,r);
u = x(1);
% y - w - z
dx(4:6) = ODE_GeneletI(x(4:6),p2,u);
dx = dx(:);
end
