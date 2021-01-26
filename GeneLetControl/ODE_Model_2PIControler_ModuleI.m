function dx = ODE_Model_2PIControler_ModuleI(t,x,p)
% Return T - RI - RA
% Parameters
p1 = p(1:6);
p2 = p(7:12);
p3 = p(13:20);
r1 = p(21); % Reference
r2 = p(22); % Reference
dx = nan(1,9);
% States
y = x(4);
% u - RA - RI
dx(1:3) = ODE_Controller(x(1:3),p1,y,r1);
% dx(4:6) = ODE_Controller(x(4:6),p1,r2,x(1));
dx(4:6) = ODE_Controller(x(4:6),p2,x(1),r1);
u = x(4);
% y - w - z
dx(7:9) = ODE_GeneletI(x(7:9),p3,u);
dx = dx(:);
end
