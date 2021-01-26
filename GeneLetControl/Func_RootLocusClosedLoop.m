function [] = Func_RootLocusClosedLoop(x,pc,ps)

[Nc,Dc] = TF_Controller(x,pc);
xc = [roots(Dc); roots(Nc)];
[Ns,Ds] = TF_Process(x,ps);
xs = [roots(Ds); roots(Ns)];
p1 = tf(conv(Nc,Ns),conv(Dc,Ds));

xcs = conv(Dc,Ds) + 1*conv(Nc,Ns);
sol = roots(xcs);

rlocus(p1)
hold on
plot(xs,xs*0,'d','Color',[160 158 158]/255)
plot(xc,xc*0,'d','Color',[248 152 56]/255,'LineWidth',2)
plot(real(sol),imag(sol),'b.')
hold off
end
function [N,D] = TF_Controller(x,p)

% Kinetic parameters
k = p(1);
th = p(2);
b = p(3);
a = p(4);
ph = p(5);
g = p(6);
ut = p(7);
% A = p(8);
I = p(9);

rA = x(1);
rI = x(2);
u = x(3);
A = x(4);

du = ut-u;

J = [-k*A/rA g*rA a*rA;
    g*rI -th*I/rI b*rI;
    a*du b*u -ut*b*rI/du];
B = [k*A;0;0];
C = [0 0 1];

[N,D] = ss2tf(J,B,C,0);

end

function [N,D] = TF_Process(x,p)

% Kinetic parameters
k = p(1);
b = p(2);
a = p(3);
ph = p(4);
g = p(5);
th = p(6);
yt = p(7);
wt = p(8);

u = x(3);
y = x(4);
w = x(5);
z = x(6);

dy = yt-y;

J = [-(a*w+b*z) a*dy b*y;
    a*w-th -(th+a*dy+g*z) g*w;
    b*z g*z -(b*y+g*w+ph)];
B = [0;0;k*u];
C = [1 0 0];

[N,D] = ss2tf(J,B,C,0);

end
