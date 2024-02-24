% Zane Grothe
% AERO 6330
% HW 4
% 3/15/22

% Function file for propagating equations of motion for trajectory
% Input:
%           t        = timespan
%           xy0      = initial conditions
%           options  = tolerances
%           mu       = nondimensional mass parameter
% Output:
%           ydot     = velocity and acceleration

function ydot = CR3BP_EOM(~,xy0,~,mu)

% Distances
d=sqrt((mu+xy0(1))^2+(xy0(2))^2);
r=sqrt((1-mu-xy0(1))^2+(xy0(2))^2);

% Outputs
ydot=[xy0(3);
      xy0(4);
      xy0(1)+2*xy0(4)+(1-mu)*(-mu-xy0(1))/(d^3)+mu*(1-mu-xy0(1))/(r^3); 
      xy0(2)-2*xy0(3)-(1-mu)*(xy0(2))/(d^3)-mu*xy0(2)/(r^3)];
end

