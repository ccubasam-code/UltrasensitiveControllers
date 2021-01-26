function dx = ODE_ControllerGene(t,x,p)
% Return T - RI - RA
% Parameters
p1 = p(1:9);
p2 = p(10:15);
r = p(16); % Reference

dx = nan(1,6);
% States
y = x(6);
% u - RA - RI
dx(1:3) = ODE_Controller(x(1:3),p1,y,r);
u = x(1);
% y - w - z
dx(4:6) = ODE_Gene(x(4:6),p2,u);
dx = dx(:);
end
