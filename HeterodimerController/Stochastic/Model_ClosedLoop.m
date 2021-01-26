function [ S, h, endSim ] = Model_ClosedLoop(Omega,p)
% Parameters
p1 = p(1:13);
p2 = p(14:17);
r = p(18); % Reference


% Rate constants 
% c = [1 d 1 d 1 d 1 d 1 d 1 d]';
c = [1 1 1 1 1 1 1 1 1 1 1 1 1 ]';
% Mass convervation conditions
% cf = @(x) (ct - x(4,:)  - (gTt-x(3,:)));
% Stoichiometry matrix

S1 = [1 -1  0   0  -1 -2  2   0  0
      0  0  1  -1  -1  0  0  -2  2
      0  0  0   0   0  1 -1   0  0
      0  0  0   0   0  0  0   1 -1
     ];

S2 = [1 -1 0  0 
      0  0 1 -1 
     ];
S = [S1 zeros(size(S1,1),size(S2,2));zeros(size(S2,1),size(S1,2)) S2];

% Rates      
h = @(x)(repmat(c,[1 size(x,2)])...
                 .*[ p1(1)*r*Omega/(r*Omega+p(12)*Omega)*ones(1,size(x,2))*Omega
                     p1(4)*x(1,:)
                     p1(2)*x(6,:)./(x(6,:) + p(13)*Omega)*Omega
                     p1(4)*x(2,:)
                     p1(3)*x(1,:).*x(2,:)/Omega
                     p1(6)*x(1,:).*(x(1,:)-1)/2.*(p1(11)*Omega-x(3,:)-x(4,:))/Omega^2
                     p1(7)*x(3,:)
                     p1(8)*x(2,:).*(x(2,:)-1)/2.*(p1(11)*Omega-x(3,:)-x(4,:))/Omega^2
                     p1(9)*x(4,:)
                     
                     p2(1)*x(3,:)
                     p2(2)*x(5,:)
                     p2(3)*x(5,:)
                     p2(4)*x(6,:)
                 
                 ]);

% Never end the simulation
endSim = @(x)(false(1,size(x,2)));


end