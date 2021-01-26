h = 3600;
uM = 10^(-6);
uMh = uM*h;
set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultAxesFontSize', 10)
% Kinetic rate of biocontroller
kc = 20.*10^(-4)*h/4;
thc = 20.*10^(-4)*h/4*1.5;
bc = 3*10^4*uMh*2;
ac = 3*10^4*uMh*2;
phc = log(2)*60/30; 
gc = 3*10^4*uMh;
ut = 0.5*2;

p1 = [kc thc bc ac phc gc ut]; % 7

% Kinetic rate of biocontroller plant
ks = kc;
bs = bc;
as = ac;
phs = phc;
gs = gc;
ths = 3.*10^(-4)*h;
yt = 0.8;
wt = yt;
r1 = 0.1;
r2 = 0.2;
r3 = 0.3;
R1 = thc*r1/kc;
R2 = thc*r2/kc;
R3 = thc*r3/kc;

p2 = [ks bs as phs gs ths yt wt]; % 8-15

p = [p1 p2 r1];
p0= p;
NN = 1000;
Color3 = [8 104 172; 67 162 202; 123 204 196]/255;

N = 7;
r_v = linspace(0,.6,N); r_v(1)=[];
for r=1:N-1
%     r
    i = 1;
    while i < 100
        nn = (rand(1,8)-0.5)*2*log10(1.2);
        p(8:15) = p0(8:15).*10.^nn;
        u_c = logspace(-6,0,NN)*p(7); %e_v(end) = 0.99999*et;
        y_c = A_F_u_I(u_c,[p(1:7) 0 r_v(r)]);

        y_s = logspace(-6,0,NN)*p(14);
        u_s = u_F_y(y_s,p(8:15));

        [u_i,y_i] = intersections(u_c,y_c,u_s,y_s,1);
        
        if (0)
            figure(2)
            plot(y_c,u_c,'Color',[248 152 56]/255,'LineWidth',2), hold on
            plot(y_s,u_s,'Color',Color3(1,:),'LineWidth',2)
            plot(y_i,u_i,'ko','LineWidth',2)
            hold off
            ylim([0 1]*ut)
            pause
        end
        if ~isempty(y_i)
            figure(7),plot(r_v(r),y_i,'k.'), hold on
            i = i+1;
        end
    end
end
hold off
