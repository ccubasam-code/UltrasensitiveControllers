function [ S, h, endSim ] = Model_Brink(Omega,p)
% Parameters
k = p(1)*Omega;
th = p(2)*Omega;
b = p(3)/Omega;
a = p(4)/Omega;
ph = p(5);
g = p(6)/Omega;
ut = p(7)*Omega;

% Rate constants 
% c = [1 d 1 d 1 d 1 d 1 d 1 d]';
c = ones(7,1);
% Mass convervation conditions
% cf = @(x) (ct - x(4,:)  - (gTt-x(3,:)));
% Stoichiometry matrix

S = [1 -1 0  0 -1 0 0
     0  0 1 -1 -1 0 0
     0  0 0  0  0 1 -1];

% Rates
              
h = @(x)(repmat(c,[1 size(x,2)])...
                 .*[ k*ones(1,size(x,2))
                     ph*x(1,:)
                     th*ones(1,size(x,2))
                     ph*x(2,:)
                     g*x(1,:).*x(2,:)
                     a*x(1,:).*(ut-x(3,:))
                     b*x(2,:).*x(3,:)
                 ]);

% Never end the simulation
endSim = @(x)(false(1,size(x,2)));


end