clear all,clc
% Exercise 4
% ----------

% This script generates a number of (Gillespie) realizations of the Goodwin
% model and plots the average and spread of the trajectories (characterized
% by the 2.5 and 97.5 percentiles) against time.
% Compare the average trajectories for p = 6 and p = 10 in the Goodwin
% model against with what you would expect from the deterministic
% simulations.


Omega = 600;

h = 3600;
uM = 10^(-6);
uMh = uM*h;
% Kinetic rate of biocontroller
kc = 21.*10^(-4)*h/4;
bc = 3*10^4*uMh/Omega;
ac = 3*10^4*uMh/Omega;
phc = log(2)*60/30; 
gc = 3*10^4*uMh/Omega;
ut = 0.5*Omega*2;

p1 = [kc bc ac phc gc ut]; % 6

% Kinetic rate of biocontroller plant
ks = kc;
bs = bc;
as = ac;
phs = phc;
gs = gc;
ths = 3.*10^(-4)*h;
yt = 0.6*Omega;
wt = yt;
r1 = 0.1*Omega;
% param = 7;
p2 = [ks bs as phs gs ths yt wt];

var = {'\kappa_c','\beta_c','\alpha_c','\phi_c','\gamma_c', 'u^{tot}', ...
    '\kappa_s','\beta_s','\alpha_s','\phi_s','\gamma_s','\theta_s','y^{tot}','w^{tot}'};

set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultAxesFontSize', 12)

p = [p1 p2 r1];

% Get reactions
[ S, h, endSim ] = Model_BioController( Omega,p );

% Initial condition
x0 = round(Omega*[0 0 0 0 0 0]');

% Simulation time
Tmax = 20;
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
tst = {'u','r_Y','r_R','y','w','z'};
figure(20)
for i=1:size(X,1)
    data = [];
    for j=1:size(X,3)
        data = [data;X(i,:,j)];
    end
    data = data'*scaleFactor;
    subplot(2,3,i)
    fanChart(Treg,data,'mean',[2.5 97.5],'alpha',.2,'colormap',{'shadesOfColor',[0 0 .8]})
    xlabel('Time (h)')
    title(tst(i))
    %ylabel('Concentration')
end