% Zane Grothe
% Auburn University
% AERO 7970
% Exam 2
% 5/5/23

% This script calculates the generalized position and momenta of a massless
% spacecraft in the circular-resticted three-body problem using Hamiltonian
% mechanics and a 2nd-order symplectic integrator (Stormer-Verlet). The
% 2nd-order integrator uses 3 steps, of which the first 2 are implicit.
% Picard iteration is used to find these values. The example below plots
% the spacecraft's trajectory in the rotating x/y plane of the Earth-Moon
% system. Several initial condition cases are given below to compare
% example trajectories. 


clear all
close all
clc


% Constants

M1=5.9722*10^24;            % Mass of primary (Earth) kg
M2=7.3420*10^22;            % Mass of secondary (Moon) kg
M=M1+M2;                    % Total mass
mu=M2/M;                    % Nondimensional mass parameter

% Initial Conditions (nondimensionalized)

%    For Earth orbit
q = [-mu;0.2;0];           % Position
p = [2.2;0.1;0];           % Momentum
%    For Lunar orbit
%q = [1-mu;0.1;0];          % Position
%p = [0.3;1;0];             % Momentum
%    For Lunar injection (almost)
%q = [-0.5;0.01;0];         % Position
%p = [0.84986;-1.34727;0];  % Momentum
%    For random fun!
%q = [rand(2,1);0];         % Position
%p = [rand(2,1);0];         % Momentum
%    Choose your own
%q = [;;0];                 % Position (ideally values between -1 and +1)
%p = [;;0];                 % Momentum

h  = 0.01;                  % Step size
qm = q;                     % q matrix to begin logging q steps
pm = p;                     % p matrix to begin logging p steps

%    Distance from Earth (nondimensionalized)
r1 = sqrt((mu+q(1))^2+q(2)^2+q(3)^2);
%    Distance from Moon (nondimensionalized)
r2 = sqrt((q(1)-(1-mu))^2+q(2)^2+q(3)^2);

%       Pseudopotential (Vbar) partial derivatives (d)
Vbard = [-q(1)+(1-mu)*(mu+q(1))/r1^3+mu*(q(1)-(1-mu))/r2^3;...
         -q(2)+(1-mu)*q(2)/r1^3+mu*q(2)/r2^3;...
         (1-mu)*q(3)/r1^3+mu*q(3)/r2^3];


%% Numerically integrate

for i = 1:1000-1        % Total steps (including initial conditions)
    
% FIRST STEP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% p(n+1/2) = p(n) - h/2 * dh/dq[ p(n+1/2) , q(n) ]
% Implicit integration involves Picard iteration
    
    tol      = 1e-9;    % Tolerance
    max_iter = 1000;    % Maximum number of Picard iterations
    iter     = 0;       % Iteration counter
    p_old    = p/2;     % Initial p(n+1/2) guess
    error    = 1;       % Initial error (large)
    
    while error > tol && iter < max_iter
              % Partial of Hamiltonian wrt q (pass in guess)
        dHdq  = -[p_old(2)-q(1)-Vbard(1);-p_old(1)-q(2)-Vbard(2);-Vbard(3)];
              % Calculate p(n+1/2)
        p_new = p - h/2 * dHdq;
              % Check error
        error = norm(abs(p_new - p_old) / abs(p_old));
              % Rename for next iteration
        p_old = p_new;
              % Next iteration
        iter  = iter + 1;
    end
    
    %  Check if iterations are converging
    if error > tol
        disp('Picard iteration has failed')
    end
    
          % Partial of Hamiltonian wrt q (pass in final step values)
    dHdq  = -[p_new(2)-q(1)-Vbard(1);-p_new(1)-q(2)-Vbard(2);-Vbard(3)];
          % Calculate final p(n+1/2) values
    pHALF = p - h/2 * dHdq;
    
    
% SECOND STEP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% q(n+1) = q(n) + h/2 * { dh/dp[ p(n+1/2) , q(n) ] + dH/dp[ p(n+1/2) ,
%          q(n+1) ] }
% Implicit integration involves Picard iteration
    
    iter     = 0;       % Iteration counter
    q_old    = q/2;     % Initial q(n+1) guess
    error    = 1;       % Initial error (large)
    
    while error > tol && iter < max_iter
                 % Partial of Hamiltonian wrt p (explicit)
        dHdp_qn  = [pHALF(1)+q(2);pHALF(2)-q(1);pHALF(3)];
                 % Partial of Hamiltonian wrt p (pass in guess)
        dHdp_qn1 = [pHALF(1)+q_old(2);pHALF(2)-q_old(1);pHALF(3)];
                 % Calculate q(n+1)
        q_new    = q + h/2 * (dHdp_qn+dHdp_qn1);
                 % Check error
        error    = norm(abs(q_new - q_old) / abs(q_old));
                 % Rename for next iteration
        q_old    = q_new;
                 % Next iteration
        iter     = iter + 1;
    end
    
    %  Check if iterations are converging
    if error > tol
        disp('Picard iteration has failed')
    end
    
             % Partial of Hamiltonian wrt p (pass in final step values)
    dHdp_qn1 = [pHALF(1)+q_new(2);pHALF(2)-q_new(1);pHALF(3)];
             % Calculate final q(n+1) values
    q        = q + h/2 * (dHdp_qn+dHdp_qn1);
             % Log new q values in q matrix
    qm       = [qm,q]; %#ok<AGROW>
    
    
% THIRD STEP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% p(n+1) = p(n+1/2) - h/2 * dh/dq[ p(n+1/2) , q(n+1) ]
% Explicit integration
    
          % Recalculate distances (pass in new q values)
    r1    = sqrt((mu+q(1))^2+q(2)^2+q(3)^2);
    r2    = sqrt((q(1)-(1-mu))^2+q(2)^2+q(3)^2);
          % Recalculate pseudopotential partial derivatives (pass in new q
          % values
    Vbard = [-q(1)+(1-mu)*(mu+q(1))/r1^3+mu*(q(1)-(1-mu))/r2^3;...
             -q(2)+(1-mu)*q(2)/r1^3+mu*q(2)/r2^3;...
             (1-mu)*q(3)/r1^3+mu*q(3)/r2^3];
          % Recalculate new partial of Hamiltonian wrt q (pass in new p and
          % q values)
    dHdq  = -[pHALF(2)-q(1)-Vbard(1);-pHALF(1)-q(2)-Vbard(2);-Vbard(3)];
          % Calculate final p(n+1) values
    p     = pHALF - h/2 * dHdq;
          % Log new p values in p matrix
    pm    = [pm,p]; %#ok<AGROW>
    
end


%% Plot
figure(1)

% Earth
plot(-mu,0,'ko','MarkerSize',8,'MarkerFaceColor','g')
hold on
% Moon
plot(1-mu,0,'ko','MarkerSize',2,'MarkerFaceColor','b')
% Trajectory
plot(qm(1,:), qm(2,:), 'k')

title('CR3BP (Earth-Moon) Trajectory (nondimensional)')
xlabel('X')
ylabel('Y')
xlim([-0.8 1.2])
ylim([-1 1])
axis square

