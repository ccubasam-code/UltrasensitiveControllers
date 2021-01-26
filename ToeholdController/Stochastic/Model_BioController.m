function [ S, h, endSim ] = Model_BioController(Omega,p)
% Parameters
pc = p(1:9);
ps = p(10:15);

ut = pc(7);
KA = pc(8);
KI = pc(9);
n = ps(6);

r = p(16); % Reference

% c = [1 d 1 d 1 d 1 d 1 d 1 d]';
c = [pc(1) pc(2) pc(6) pc(5) pc(5) pc(4) pc(3) ...
    ps(4) ps(5) ps(1) ps(3) ps(4) ps(5)]';
% Mass convervation conditions
% cf = @(x) (ct - x(4,:)  - (gTt-x(3,:)));
% Stoichiometry matrix
A = zeros(3,7);
A(1,[1 3 4 6]) = [1 -1 -1 -1];
A(2,[2 3 5 7]) = [1 -1 -1 -1]; 
A(3,[6 7]) = [1 -1]; 

B = zeros(3,6);
B(1,[1 2]) = [1 -1];
B(2,[3 4]) = [1 -1];
B(3,[5 6]) = [1 -1];

S = [A B*0;0*A B];
% Rates
              
h = @(x)(repmat(c,[1 size(x,2)])...
                 .*[x(6,:)./(x(6,:)+KA)
                    ones(1,size(x,2))*r/(r+KI)
                    x(1,:).*x(2,:)
                    x(1,:)
                    x(2,:)
                    x(1,:).*(ut-x(3,:))
                    x(2,:).*x(3,:)
                    x(3,:)
                    x(4,:)
                    KA^n./(x(4,:).^n+KA^n)
                    x(5,:)
                    x(5,:)
                    x(6,:)
                 ]);

% Never end the simulation
endSim = @(x)(false(1,size(x,2)));


end