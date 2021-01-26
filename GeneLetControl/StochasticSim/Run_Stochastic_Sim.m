function [X,T,x,t] = Run_Stochastic_Sim(Omega,p,x,t,Tmax)

X = nan(size(x,1),100000);
T = nan(1,100000);
% Write the initial state to the first "record"
X(:,1) = x;
T(1) = t;

[ S, h, endSim ] = Model_BioController( Omega,p );
step = 1;
while t<Tmax
    
    hs = h(x);
    hscs = cumsum(hs);
    rateall = sum(hs);
    tau = -log(rand)/rateall;
    r = rand*sum(hs);
    j = find(r<=hscs,1,'first');
    x = x + S(:,j);
    t = t+tau;
    step = step + 1;
    X(:,step) = x;
    T(step) = t;
end