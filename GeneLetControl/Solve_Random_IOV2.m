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
ut = 0.5;
j = 1;
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
NN = 2000;
Color3 = [8 104 172; 67 162 202; 123 204 196]/255;
N = 1000;

hFig=figure(1);
set(hFig,'Units','inches', 'Position', [0 6.5 4.5 1.5])
utv = [0.5 1 2]*ut;
% utv = 0.5*ut;
for j = 3
    p1 = [kc thc bc ac phc gc utv(j)]; % 7
    p = [p1 p2 r1];
Data = nan(N,2);
Datap = nan(N,18);
i = 1;
    while i < N+1
        r_v = rand(1)*0.6;
        nn = (rand(1,8)-0.5)*2*log10(1.5);
        p(8:15) = p0(8:15).*10.^nn;
        u_c = logspace(-6,0,NN)*p(7); %e_v(end) = 0.99999*et;
        y_c = A_F_u_I(u_c,[p(1:7) 0 r_v]);

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
%             figure(7),plot(r_v,y_i,'k.'), hold on
            Data(i,:) =[ thc*r_v/kc y_i];
            Datap(i,:) =[p r_v y_i];
            i = i+1;
        end
    end

    subplot(1,3,j),plot(Data(:,1),Data(:,2),'k.','LineWidth',.1)
    xlim([0 0.8])
    ylim([0 0.8])
    xlabel('r (\mu M)')
%     ylabel('y (\mu M)')
    ax = gca;
    ax.XTick = [0  0.8];
    ax.YTick = [0 0.8];
end


%%
idx1 = find(Data(:,2)./Data(:,1)>1);
idx2 = find(Data(:,2)./Data(:,1)<0.8);
figure(4)
subplot(1,4,2)
    
plot(Data(:,2)./Data(:,1))

ii = idx2(2);
jj = idx1(2);
subplot(1,4,1)
plot(Data(:,1),Data(:,2),'k.'), hold on
plot(Data(ii,1),Data(ii,2),'ro')
plot(Data(jj,1),Data(jj,2),'bo'),hold off
xlim([0 0.8])
ylim([0 0.8])

if (0)
    figure(1),subplot(1,3,j)
    plot(Data(:,1),Data(:,2),'k.'), hold on
% plot(Data(ii,1),Data(ii,2),'go')
plot(Data(jj,1),Data(jj,2),'bo'),hold off
xlim([0 0.8])
    ylim([0 0.8])
    xlabel('r (\mu M)')
%     ylabel('y (\mu M)')
    ax = gca;
    ax.XTick = [0  0.8];
    ax.YTick = [0 0.8];
end

figure(4)
p = Datap(ii,[1:15,17]);
% p(7)=2*p(7);
NN = 10000;
u_c = logspace(-6,0,NN)*p(7); %e_v(end) = 0.99999*et;
y_c = A_F_u_I(u_c,[p(1:7) 0 p(end)]);

y_s = logspace(-6,0,NN)*p(14);
u_s = u_F_y(y_s,p(8:15));

[u_i,y_i] = intersections(u_c,y_c,u_s,y_s,1);

subplot(1,4,3)
plot(y_c,u_c,'Color',[248 152 56]/255,'LineWidth',2), hold on
plot(y_s,u_s,'Color',Color3(1,:),'LineWidth',2)
plot(y_i,u_i,'ro','LineWidth',2)
hold off
xlim([0 1])
ylim([0 1]*p(7))
axis square

p = Datap(jj,[1:15,17]);
% p(7)=2*p(7);
u_c = logspace(-6,0,NN)*p(7); %e_v(end) = 0.99999*et;
y_c = A_F_u_I(u_c,[p(1:7) 0 p(end)]);

y_s = logspace(-6,0,NN)*p(14);
u_s = u_F_y(y_s,p(8:15));

[u_i,y_i] = intersections(u_c,y_c,u_s,y_s,1);

subplot(1,4,4)
plot(y_c,u_c,'Color',[248 152 56]/255,'LineWidth',2), hold on
plot(y_s,u_s,'Color',Color3(1,:),'LineWidth',2)
plot(y_i,u_i,'bo','LineWidth',2)
hold off
xlim([0 1])
ylim([0 1]*p(7))
axis square
            