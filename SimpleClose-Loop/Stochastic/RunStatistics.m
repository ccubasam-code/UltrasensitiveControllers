% Solve RNA clock
%close all,clc, clear all
% Bistable parametrs
h = 3600;
uM = 10^(-6);
uMh = uM*h;
F = 1;
% Kinetic rate of biocontroller
kc = 20.*10^(-4)*h/4*F;
thc = 20.*10^(-4)*h/4*1.5*F;
bc = 3*10^(4)*uMh*2*1;
ac = 3*10^(4)*uMh*2*1;
phc = log(2)*60/30; 
gc = 3*10^4*uMh;
ut = 0.5*1;

pI = [kc thc bc ac phc gc ut]; % 5-11

% Kinetic rate of biocontroller plant
as = 20.*10^(-4)*h*.5/1;
ks = .2/1; % /2
phs = log(2)*60/3/1;
trs = .23*60/4*1; %*2
ds = log(2)*60/30;
ns = 2;
yt = 0.6;
pS = [as ks phs trs ds ns]; % 12-17

r1 = 0.1/1.5*2; % 2 3 4
R1 = thc*r1/kc;



c1 = 1;% proportional 
c2 = 1;% Integral
p = [pI pS r1];

Omega = 600/1;

set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultAxesFontSize', 12)


% Get reactions
[ S, h, endSim ] = Model_Feedback( Omega,p );

% Initial condition
x0 = round(Omega*[0. 0.  0 0 0 ]')*0;

% Simulation time
Tmax = 20;
%///////////////////////////////

set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultAxesFontSize', 12)

%% Set up simulation (multiple trajectories)

reset(RandStream.getGlobalStream);

% Number of trajectories to simulate
numReals = 500; % 1000

% Initial state and time
t = repmat(0,[1 numReals]);
x = repmat(x0,[1 numReals]);

% Number of steps done
idx = 0;

% Set up arrays to record trajectories
tsampleIdx = ones(1,numReals);
deltaSample = 0.1;
X = nan(size(x,1),ceil(Tmax/deltaSample+1),size(x,2));
Treg = (0:size(X,2)-1)*deltaSample;
T = nan(1,size(X,2),size(x,2));
X(:,1,:) = x;
T(1,1,:) = t;

%% Actual simulation

while any(t<Tmax)
    updateThese = t<Tmax;
    
    if mod(idx,100)==0
        disp([num2str(sum(updateThese)) ' to update, median(t) = ' num2str(median(t(updateThese)))]);
    end
    idx = idx + 1;
    
    [deltaT, deltaX] = stepGillespie(x(:,updateThese),h,S);
    
    lastx = x;
    lastt = t;

    % Update state and time
    x(:,updateThese) = x(:,updateThese) + deltaX;
    t(updateThese) = t(updateThese) + deltaT;
    
    
    while(true)
        saveThese = t>(1+eps)*tsampleIdx*deltaSample & tsampleIdx+1<=size(X,2);
        if(~any(saveThese))
            break;
        end
        tsampleIdx(saveThese) = tsampleIdx(saveThese)+1;
        for saveThis = find(saveThese)
            X(:,tsampleIdx(saveThis),saveThis) = lastx(:,saveThis);
            T(1,tsampleIdx(saveThis),saveThis) = lastt(saveThis);
        end
    end
    
    % End some simulations
    endThese = endSim(x);
    x(:,endThese) = nan;
    t(endThese)=nan;
    
end

%figure;

% Molecule Numbers
%plot(Treg,squeeze(X(1,:,:)));

% Concentrations
%plot(Treg,squeeze(X(1,:,:))/Omega);



%%
scaleFactor = 1/Omega;
tst = {'u','r_y','r_r','r_m','y'};
figure(20)
nf = [1 2 3 4 5 6 7];
for i=1:size(X,1)
    data = [];
    for j=1:size(X,3)
        data = [data;X(i,:,j)];
    end
    data = data'*scaleFactor;
    subplot(3,3,nf(i))
    fanChart(Treg,data,'mean',[2.5 97.5],'alpha',.2,'colormap',{'shadesOfColor',[0 0 .8]})
    xlabel('Time (h)')
    title(tst(i))
    %ylabel('Concentration')
    hold off
end

%%
set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultAxesFontSize', 12)
tst = {'u','r_y','r_r','r_m','y'};
hFig=figure(22);
set(hFig,'Units','inches', 'Position', [0 100 7/2 2])

i1 = 1;
i2 = 5;
data1 = [];data2 = [];
for j=1:size(X,3)
    data1 = [data1;X(i1,:,j)];
    data2 = [data2;X(i2,:,j)];
end
data1 = data1'*scaleFactor;
% data2 = data2';
data2 = data2'*scaleFactor;
% subplot(1,2,1)
% fanChart(Treg,data1,'mean',[2.5 97.5],'alpha',.2,'colormap',{'shadesOfColor',[0 0 .8]})
% xlabel('Time (h)')
% % title(tst(i))
% ylabel('\mu M')
% hold off
% ax = gca;
% ax.XTick = [0  r1 yt];
% ax.YTick = [0 0.2 .4];
% hold on
subplot(1,2,2)
fanChart(Treg,data2,'mean',[2.5 97.5],'alpha',.2,'colormap',{'shadesOfColor',[0.8 0. 0]})

xlabel('Time (h)')
% title(tst(i))
ylabel('\mu M')
hold off
ax = gca;
% ax.XTick = [0  r1 yt];
% ax.YTick = [0 0.2 .4];

