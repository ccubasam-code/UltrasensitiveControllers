h = 3600;
uM = 10^(-6);
uMh = uM*h;
set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultAxesFontSize', 10)
% Kinetic rate of biocontroller
kc = 20.*10^(-4)*h/4;
thc = 20.*10^(-4)*h/4;
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

ff = 1.5;

p2 = [ks bs as phs gs ths yt wt]; % 8-15

p = [p1 p2 r1];
p0= p;
NN = 2000;
Color3 = [160 158 158; 67 162 202; 123 204 196]/255;
N = 10000;


% utv = [0.5 1 2]*ut;
utv = ut;
for j = 1:1
    p1 = [kc thc bc ac phc gc utv(j)]; % 7
    p = [p1 p2 r1];
    
    DataO = nan(N,2);
    i = 1;
    while i < N+1
        nn = (rand(1,8)-0.5)*2*log10(ff);
        p(8:15) = p0(8:15).*10.^nn;

        y_s = rand(1)*p(14);
        u_s = u_F_y(y_s,p(8:15));
%         if imag(u_s) ~=0
            DataO(i,:) =[ 1/(1+real(u_s)/0.2) y_s];
            i = i+1;
%         end
        
    end
%

Data = nan(N,3);
Datap = nan(N,18);
i = 1;
    while i < N+1
        tmp = 0;
        r_v = rand(1)*1.5;
        p(16) = r_v;
        nn = (rand(1,15)-0.5)*2*log10(ff); % 1.2
        p(1:15) = p0(1:15).*10.^nn;
        u_c = logspace(-6,0,NN)*p(7); %e_v(end) = 0.99999*et;
        y_c = A_F_u_I(u_c,[p(1:7) 0 r_v]);

        y_s = linspace(0,1,NN)*p(14);
        u_s = u_F_y(y_s,p(8:15));
        idx = find (y_c<0);u_c(idx) = [];y_c(idx)=[];
        
        [u_i,y_i] = intersections(u_c,y_c,u_s,y_s,1);
        
             
        if ~isempty(y_i)
            
            [~,eigJ,~] = ComputeStability(p,u_i,y_i);
            if sum(real(eigJ)>0)>1 % Red
                tmp = 1;
            end
        
%             figure(7),plot(r_v,y_i,'k.'), hold on
            Data(i,:) =[ thc*r_v/kc y_i tmp];
            Datap(i,:) =[p r_v y_i];      
            i = i+1;
           
        end
    end
% end
%
% Color = [153 153 153;155 166 182; 5 32 73]/255;

Color = [163 218 222;5 32 73; 200 200 200]/255;

hFig=figure(1);
set(hFig,'Units','inches', 'Position', [0 6.5 4.5 1.5])

% hFig=figure(2);
% set(hFig,'Units','inches', 'Position', [0 6.5 4.5 1.5])


XX = 0.1;

idx0 = find(Data(:,3)==0); % idx0
idx1 = find(Data(:,3)==1); % idx1

[size(idx0) size(idx1)]

figure(1)
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

nbins = [0:0.05:2];


hFig=figure(10);
set(hFig,'Units','inches', 'Position', [0 6.5 4.5 1.5])

subplot(1,3,2)
[h2] = histogram(real(Data(:,2))./real(Data(:,1)),nbins,'Normalization','probability');
hold on;
[h1] = histogram(real(DataO(:,2))./real(DataO(:,1)),nbins,'Normalization','probability');
h1.EdgeColor = Color(1,:);
h2.EdgeColor = Color(2,:);
h1.FaceColor = Color(1,:);
h2.FaceColor = Color(2,:);
hold off
xlabel(' ')


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

size(idx1,1)
end
%%
idx1 = find(Data(:,2)./Data(:,1)>1);
idx2 = find((Data(:,2)./Data(:,1)>.9) & Data(:,1)>0.5);
figure(4)
subplot(1,5,2)
    
plot(Data(:,2)./Data(:,1))

subplot(1,5,5)
plot(Data(:,3),'k.')
ylim([0 1])

ii = idx1(1);
jj = idx2(1);
subplot(1,5,1)
plot(Data(:,1),Data(:,2),'k.'), hold on
plot(Data(ii,1),Data(ii,2),'ro')
plot(Data(jj,1),Data(jj,2),'bo')
hold off
xlim([0 0.8])
ylim([0 0.8])

if (0)
    figure(3),subplot(1,3,j)
   subplot(1,3,j),pp1=plot( DataO(:,1),DataO(:,2),'.','Color',Color(1,:)), hold on
   pp1.Color(4)=XX;
   subplot(1,3,j),pp2=plot(Data(:,1),Data(:,2),'.','Color',Color(2,:))
   pp2.Color(4)=XX;
   plot(Data(jj,1),Data(jj,2),'ro'),hold off
   xlim([0 0.8])
    ylim([0 0.8])
    xlabel('r (\mu M)')
%     ylabel('y (\mu M)')
    ax = gca;
    ax.XTick = [0  0.8];
    ax.YTick = [0 0.8];
end

hFig=figure(4);
set(hFig,'Units','inches', 'Position', [0 3.5 6 1.5])
p = Datap(ii,[1:15,17]);
% p(7)=2*p(7);
NN = 1000;
u_c = logspace(-6,0,NN)*p(7); %e_v(end) = 0.99999*et;
y_c = A_F_u_I(u_c,[p(1:7) 0 p(end)]);

y_s = logspace(-6,0,NN)*p(14);
u_s = u_F_y(y_s,p(8:15));

[u_i,y_i] = intersections(u_c,y_c,u_s,y_s,1);

subplot(1,5,3)
plot(y_c,u_c,'Color',[248 152 56]/255,'LineWidth',2), hold on
plot(y_s,u_s,'Color',[160 158 158]/255,'LineWidth',2)
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

subplot(1,5,4)
plot(y_c,u_c,'Color',[248 152 56]/255,'LineWidth',2), hold on
plot(y_s,u_s,'Color',[160 158 158]/255,'LineWidth',2)
plot(y_i,u_i,'bo','LineWidth',2)
hold off
xlim([0 1])
ylim([0 1]*p(7))
axis square
            