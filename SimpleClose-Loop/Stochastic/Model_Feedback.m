function [ S, h, endSim ] = Model_Feedback(Omega,p)
% Parameters
k = p(1);
th = p(2);
b = p(3)/Omega;
a = p(4)/Omega;
ph = p(5);
g = p(6)/Omega;
ut = p(7)*Omega;

as = p(8)*Omega;
kus = p(9)*Omega;
phs = p(10);
trs = p(11);
ds = p(12);
ns = p(13);
r = p(14)*Omega;
% c = [1 d 1 d 1 d 1 d 1 d 1 d]';
c = ones(1,11)';
% Mass convervation conditions
% cf = @(x) (ct - x(4,:)  - (gTt-x(3,:)));
% Stoichiometry matrix
S1 = [1  -1  0  0  0  0  0;
      -1  0  1 -1  0  0 -1;
      0  -1  0  0  1 -1 -1];
S2 = [1 -1 0  0;
      0  0 1 -1];

S = [S1 zeros(size(S1,1),size(S2,2));zeros(size(S2,1),size(S1,2)) S2];

% Rates
              
h = @(x)(repmat(c,[1 size(x,2)])...
                 .*[ a*(ut-x(1,:)).*x(2,:)
                     b*x(1,:).*x(3,:)
                     k*x(5,:)
                     ph*x(2,:)
                     th*r*ones(1,size(x,2))
                     ph*x(3,:)
                     g*x(2,:).*x(3,:)
                     as*kus^ns./(x(1,:).^ns+kus^ns);
                     phs*x(4,:)
                     trs*x(4,:)
                     ds*x(5,:)
                 ]);

% Never end the simulation
endSim = @(x)(false(1,size(x,2)));


end