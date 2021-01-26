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
A = 0.1;
I = 0.1;

kc = 20.*10^(-4)*h/4*A;
thc = 20.*10^(-4)*h/4*1.1*I;
bc = 3*10^4*uMh*2;
ac = 3*10^4*uMh*2;
phc = log(2)*60/30; % 3 minutes 
gc = 3*10^4*uMh;
ut = 0.5;

nnv = 6;
Omega = 600/nnv;
p = [kc thc bc ac phc gc ut]; % 6
%///////////////////////////////


% Get reactions
[ S, h, endSim ] = Model_Brink( Omega,p );

% Initial condition
x0 = round(Omega*[0. .0 0.0]');

% Simulation time
Tmax = 5*4;
%///////////////////////////////

i1 = 1;
i2 = 3;
N = 500;

[t, x, z] = ComputeStatistics(p,x0,Tmax,N,Omega,i1,i2);

%%
zm = 1;
xm = 0.6;
figure(1)
subplot(2,2,1),plot(t,x,'b')
ylim([0 xm])

nbins = [1:6:Omega]/Omega;
[h1 bins1]= hist(x(end,:),nbins);
subplot(2,2,2),barh(bins1,h1/N);
ylim([0 xm])


subplot(2,2,3),plot(t,z,'b')
ylim([0 zm])

nbins = [1:20:Omega]/Omega;
[h2 bins2]= hist(z(end,:),nbins);
subplot(2,2,4),barh(bins2,h2/N);
ylim([0 zm])

figure(2)

subplot(2,2,1),fanChart(t,x,'mean',[2.5 97.5],'alpha',.2,'colormap',{'shadesOfColor',[0 0 .8]})
ylim([0 xm])
subplot(2,2,2),barh(bins1,h1/N);
ylim([0 xm])

subplot(2,2,3),fanChart(t,z,'mean',[2.5 97.5],'alpha',.2,'colormap',{'shadesOfColor',[0.8 0. 0]})
ylim([0 zm])
subplot(2,2,4),barh(bins2,h2/N);
ylim([0 zm])

%% Sensitivity analysis

vec = [0.5 0.8 0.9 1 1.1 1.2 2]*1.1;
% vec = [0.5 1 2];
p0 = p;
i1 = 3;
i2 = 3;
N = 500*4;
nbins = [1:round(10/nnv):Omega]/Omega;

hFig=figure(101);
set(hFig,'Units','inches', 'Position', [0 1 7/2 3*2])

hFig=figure(102);
set(hFig,'Units','inches', 'Position', [0 10 7 2])


l = size(vec,2);
zm1 = 0.6;
zm2 = 0.6;
for i=1:l
    p0(1) = p(1)*vec(i);
    [t, x, z] = ComputeStatistics(p0,x0,Tmax,N,Omega,i1,i2);
    [h bins]= hist(z(end,:),nbins);
    figure(101),subplot(l,1,i)
    bar(bins,h/N);
    xlim([0 zm1])
    ylim([0 zm2])
    
    figure(102),subplot(1,l,i)
    barh(bins,h/N);
    ylim([0 zm1])
    xlim([0 zm2])
end