function [ S, h, endSim ] = Model_BioController(Omega,p)
% Parameters
p1 = p(1:7);
p2 = p(8:15);

ut = p(7);
yt = p(14);
wt = p(15);
r = p(16); % Reference

% c = [1 d 1 d 1 d 1 d 1 d 1 d]';
c = [p1(1) p1(2) p1(5) p1(4) p1(4) p1(3) p1(2) ...
    p2(1) p2(5) p2(4) p2(6) p2(3) p2(2)]';
% Mass convervation conditions
% cf = @(x) (ct - x(4,:)  - (gTt-x(3,:)));
% Stoichiometry matrix
A = zeros(3,7);
A(1,6) = 1; A(1,7) = -1;
A(2,1) = 1; A(2,3) = -1;A(2,4) = -1;A(2,6) = -1;
A(3,2) = 1; A(3,3) = -1; A(3,5) = -1;A(3,7) = -1;

B = zeros(3,6);
B(1,5) = 1; B(1,6) = -1;
B(2,2) = -1;B(2,4) = 1;B(2,5) = -1;
B(3,1) = 1;B(3,2) = -1;B(3,3) = -1;B(3,6) = -1;

S = [A B*0;0*A B];
% Rates
              
h = @(x)(repmat(c,[1 size(x,2)])...
                 .*[x(4,:)
                    ones(1,size(x,2))*r
                    x(2,:).*x(3,:)
                    x(2,:)
                    x(3,:)
                    x(2,:).*(ut-x(1,:))
                    x(3,:).*x(1,:)
                    x(1,:)
                    x(5,:).*x(6,:)
                    x(6,:)
                    (wt-x(5,:)-x(4,:))
                    x(5,:).*(yt-x(4,:))
                    x(4,:).*x(6,:)
                 ]);

% Never end the simulation
endSim = @(x)(false(1,size(x,2)));


end