function dx = ODE_ControllerGene(t,x,p)
% Return T - RI - RA
% Parameters
p1 = p(1:12);
p2 = p(13:16);
r = p(17); % Reference

dx = nan(1,6);
% States
y = x(6);
% u - RA - RI
dx(1:4) = ODE_Controller(x(1:4),p1,r,y);
uA = x(3);
% y - w - z
dx(5:6) = ODE_Gene(x(5:6),p2,uA);
dx = dx(:);
end
