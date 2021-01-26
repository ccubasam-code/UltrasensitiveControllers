function [Treg, data1, data2] = ComputeStatistics(p,x0,Tmax,N,Omega,i1,i2)

% p: parameters of the model
% x0: initial condition, normally all zeros
% Tmax: total time of simulations
% N: Sample size
% i1: state to the whole system

% It returns the stochatic simulation of the state i1 for N sample

% Get reactions
[ S, h, endSim ] = Model_Brink( Omega,p );

%///////////////////////////////


%% Set up simulation (multiple trajectories)

reset(RandStream.getGlobalStream);

% Number of trajectories to simulate
numReals = N;

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
    
%     if mod(idx,100)==0
%         disp([num2str(sum(updateThese)) ' to update, median(t) = ' num2str(median(t(updateThese)))]);
%     end
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

%% Convert data into a column matrices
data1 = [];
data2 = [];
for j=1:size(X,3)
    data1 = [data1;X(i1,:,j)];
    data2 = [data2;X(i2,:,j)];
end
data1 = data1'/Omega;
data2 = data2'/Omega;



