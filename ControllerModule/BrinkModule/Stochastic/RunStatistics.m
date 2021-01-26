clear all,clc
close all
% Exercise 4
% ----------

% This script generates a number of (Gillespie) realizations of the Goodwin
% model and plots the average and spread of the trajectories (characterized
% by the 2.5 and 97.5 percentiles) against time.
% Compare the average trajectories for p = 6 and p = 10 in the Goodwin
% model against with what you would expect from the deterministic
% simulations.
h = 3600;
uM = 10^(-6);
uMh = uM*h;
set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultAxesFontSize', 10)
% Kinetic rate of biocontroller
A = 0.15;
I = 0.1;

kc = 20.*10^(-4)*h/4*A;
thc = 20.*10^(-4)*h/4*1.1*I; % 1.5
bc = 3*10^4*uMh*2;
ac = 3*10^4*uMh*2;
phc = log(2)*60/30; % 3 minutes 
gc = 3*10^4*uMh;
ut = 0.5;


Omega = 600;
p = [kc thc bc ac phc gc ut]; % 6
%///////////////////////////////


% Get reactions
[ S, h, endSim ] = Model_Brink( Omega,p );

% Initial condition
x0 = round(Omega*[0. .0 0.0]');

% Simulation time
Tmax = 5*4;
%///////////////////////////////


%% Set up simulation (multiple trajectories)

reset(RandStream.getGlobalStream);

% Number of trajectories to simulate
numReals = 500;

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
%     pause
    
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



%% Plot statistics
%figure(1);
colors = {'b','g','r'};

scaleFactor = 1;

% Comment out for molecule numbers
scaleFactor = 1/Omega;
if (0)
for idx=1:size(X,1)
    trace = nanmean(X(idx,:,:),3);
    lower = prctile(X(idx,:,:),2.5,3);
    upper = prctile(X(idx,:,:),97.5,3);
    
    
    hp=plot(Treg,trace*scaleFactor,'LineWidth',2,'Color',colors{idx});
    hold off;
    hold on;
    plot(Treg,[lower; upper]*scaleFactor,'--','Color',get(hp,'Color'));
    hold all;
end
hold off;
end
%%
if(1)
tst = {'r_A','r_I','u'};
figure(20)
for i=1:size(X,1)
    data1 = [];
    for j=1:size(X,3)
        data1 = [data1;X(i,:,j)];
    end
    data1 = data1'*scaleFactor;
    subplot(2,3,i)
    fanChart(Treg,data1,'mean',[2.5 97.5],'alpha',.2,'colormap',{'shadesOfColor',[0 0 .8]})
    xlabel('Time (h)')
    title(tst(i))
    %ylabel('Concentration')
end
hold off
end
%%
%%
set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultAxesFontSize', 12)
tst = {'x_1','y_1','w_1','x_2','y_2','w_2','u'};
hFig=figure(22);
set(hFig,'Units','inches', 'Position', [0 100 6 3]*1.3)

i1 = 1;
i2 = 3;
data1 = [];data2 = [];
for j=1:size(X,3)
    data1 = [data1;X(i1,:,j)];
    data2 = [data2;X(i2,:,j)];
end
data1 = data1'*scaleFactor;
data2 = data2'*scaleFactor;
subplot(2,1,1),fanChart(Treg,data1,'mean',[2.5 97.5],'alpha',.2,'colormap',{'shadesOfColor',[0 0 .8]})
% hold on
subplot(2,1,2),fanChart(Treg,data2,'mean',[2.5 97.5],'alpha',.2,'colormap',{'shadesOfColor',[0.8 0. 0]})

xlabel('Time (h)')
% title(tst(i))
ylabel('\mu M')
% hold off
ax = gca;
% ax.XTick = [0  r1 yt];
% ax.YTick = [0 0.5 1 1.5];
%%
% figure(22),print -depsc2 StochSimRepre.eps
% figure(2),print -depsc2 StochSimRepreSingle.eps

