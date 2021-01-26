wh = 3600;
uM = 10^(-6);
uMh = uM*h;
% Kinetic rate of biocontroller
kc = 20.*10^(-4)*h/4/4;
thc = 20.*10^(-4)*h/4/4;
bc = 3*10^4*uMh*2*10;
ac = 3*10^4*uMh*2*10;
phc = log(2)*60/12; % 3 minutes
dc = log(2)*60/12; % 3 minutes
gc = 3*10^4*uMh/1;
ut = 0.5*1;
KAc = 0.4;
KIc = 0.3;

r1 = 0.09;
j = 1;

p1 = [kc thc bc ac phc gc ut KAc KIc ]; % 1-9

% Kinetic rate of biocontroller plant
as = 20.*10^(-4)*h*.5;
ks = .2;
phs = log(2)*60/3;
trs = .23*60/4;
ds = log(2)*60/30;
ns = 2;

p2 = [as ks phs trs ds ns ]; % 10-15

p = [p1 p2 r1];
p0= p;
NN = 5000;
Color3 = [160 158 158; 67 162 202; 123 204 196]/255;
N = 10000;



DataO = nan(N,2);
i = 1;
    while i < N+1
        y_s = rand(1)*0.8;
        nn = (rand(1,6)-0.5)*2*log10(1.2); % 1.2
        p(10:15) = p0(10:15).*10.^nn;

        u_s = u_F_y(y_s,p(10:15));
%         if imag(u_s) ~=0
            DataO(i,:) =[ 1/(1+real(u_s)/0.1) y_s];
            i = i+1;
%         end
        
    end
%
Data = nan(N,2);
Datap = nan(N,18);
i = 1;
    while i < N+1
        r_v = rand(1)*0.6;
        nn = (rand(1,6)-0.5)*2*log10(1.2); % 1.2
        p(10:15) = p0(10:15).*10.^nn;
        
        u_c = linspace(0,1,NN)*p(7); %e_v(end) = 0.99999*et;
        y_c = A_F_u_I(u_c,[p(1:9) 0 r_v]);
        idx = find (y_c<0);u_c(idx) = [];y_c(idx)=[];

        y_s = logspace(-6,0,NN/50);
        u_s = u_F_y(y_s,p(10:15));

        [u_i,y_i] = intersections(u_c,y_c,u_s,y_s,1);
        
        if (0)
            figure(10)
            plot(y_c,u_c,'Color',[248 152 56]/255,'LineWidth',2), hold on
            plot(y_s,u_s,'Color',Color3(1,:),'LineWidth',2)
            plot(y_i,u_i,'ko','LineWidth',2)
            hold off
            ylim([0 1]*ut)
            pause
        end
        if length(y_i)==1 
%             figure(7),plot(r_v,y_i,'k.'), hold on
            tmp = KAc/(kc/thc*(1+KIc/real(r_v))-1);
%             if imag(tmp) ~=0 
                Data(i,:) =[ tmp y_i];
                Datap(i,:) =[p r_v y_i];
                i = i+1;
%             end
        end
    end

%%
Color = [163 218 222;5 32 73]/255;
hFig=figure(1);
set(hFig,'Units','inches', 'Position', [0 6.5 4.5 1.5])

XX = 0.1;
subplot(1,3,j),pp1=plot( DataO(:,1),DataO(:,2),'.','Color',Color(1,:)), hold on
pp1.Color(4)=XX;
subplot(1,3,j),pp2=plot(Data(:,1),Data(:,2),'.','Color',Color(2,:))
pp2.Color(4)=XX;
xlim([0 0.8])
ylim([0 0.8])

xlabel('r (\mu M)')
ax = gca;
ax.XTick = [0  0.8];
ax.YTick = [0 0.8];
hold off

%
nbins = [0:0.05:2];

figure(100)
[h2] = histogram(real(Data(:,2))./real(Data(:,1)),nbins,'Normalization','probability');
hold on;
[h1] = histogram(real(DataO(:,2))./real(DataO(:,1)),nbins,'Normalization','probability');
h1.EdgeColor = Color(1,:);
h2.EdgeColor = Color(2,:);
h1.FaceColor = Color(1,:);
h2.FaceColor = Color(2,:);

figure(1)
subplot(1,3,2)
plot(h1.BinEdges(1:end-1)+0.05/2, h1.Values, 'Color',Color(1,:),'LineWidth',2);
hold on;
plot(h2.BinEdges(1:end-1)+0.05/2, h2.Values, 'Color',Color(2,:),'LineWidth',2);
% grid on;
hold off

hold off
xlim([0. 2])
ylim([0 1])
xlabel('a')
ax = gca;
ax.YTick = [0 1];
figure(2)
plot((DataO(:,2))./(DataO(:,1)))


idx2 = find(Data(:,2)./Data(:,1)<0.8);
figure(4)
subplot(1,4,2)
    
plot(real(Data(:,2))./real(Data(:,1)))

ii = idx2(1);


figure(1)
subplot(1,4,1)
XX = 0.1;
subplot(1,3,j),pp1=plot( DataO(:,1),DataO(:,2),'.','Color',Color(1,:)), hold on
pp1.Color(4)=XX;
subplot(1,3,j),pp2=plot(Data(:,1),Data(:,2),'.','Color',Color(2,:))
pp2.Color(4)=XX;

plot(Data(ii,1),Data(ii,2),'ro')
hold off

xlim([0 0.8])
ylim([0 0.8])

xlabel('r (\mu M)')
ax = gca;
ax.XTick = [0  0.8];
ax.YTick = [0 0.8];
hold off


hFig=figure(5);
set(hFig,'Units','inches', 'Position', [0 2.5 4.5 1.5])
p = Datap(ii,[1:15,17]);
% p(7)=2*p(7);
NN = 5000;
u_c = logspace(-6,0,NN)*p(7); %e_v(end) = 0.99999*et;
y_c = A_F_u_I(u_c,[p(1:9) 0 p(end)]);
idx = find (y_c<0);u_c(idx) = [];y_c(idx)=[];

y_s = logspace(-6,0,NN/50);
u_s = u_F_y(y_s,p(10:15));

[u_i,y_i] = intersections(u_c,y_c,u_s,y_s,1);

subplot(1,4,3)
plot(y_c,u_c,'Color',[248 152 56]/255,'LineWidth',2), hold on
plot(y_s,u_s,'Color',[160 158 158]/255,'LineWidth',2)
plot(y_i,u_i,'ro','LineWidth',2)
hold off
xlim([0 1])
ylim([0 1]*p(7))
axis square

p = Datap(ii,[1:15,17]);
% p(7)=2*p(7);
u_c = logspace(-6,0,NN)*p(7); %e_v(end) = 0.99999*et;
y_c = A_F_u_I(u_c,[p(1:9) 0 p(end)]);
idx = find (y_c<0);u_c(idx) = [];y_c(idx)=[];

y_s = logspace(-6,0,NN/50);
u_s = u_F_y(y_s,p(10:15));

[u_i,y_i] = intersections(u_c,y_c,u_s,y_s,1);

subplot(1,4,4)
plot(y_c,u_c,'Color',[248 152 56]/255,'LineWidth',2), hold on
plot(y_s,u_s,'Color',[160 158 158]/255,'LineWidth',2)
plot(y_i,u_i,'bo','LineWidth',2)
hold off
xlim([0 1])
ylim([0 1]*p(7))
axis square


            