%% PROBLEM 3 - Arithmetic Solution
close all;
clear;
% *** Given Problem Parameters ***
Qf = [100 0; 0 100];
Q = eye(2); R = 1;
T = 100; w = 2*pi;

repelemFactor = 2;
line_width = 1;
x10 = 50; x20 = 140;

%% Solution to time variable system
disp('************************************************************');
disp('Working on time varying problem...');
% ***State equations***
syms t;

A = [0, cos(w*t); sin(w*t), 0];
AT = [0, sin(w*t); cos(w*t), 0];
B = [0; 1];

% *** Ricatti Costate solutions ***
syms p1(t) p2(t) p(t);
P = [p1(t), p(t); p(t), p2(t)];
DP = -(P*A+AT*P- P*(B*B')*P+Q);

eqP1 = diff(p1(t),t) == DP(1,1);
eqP2 = diff(p2(t),t) == DP(2,2);
eqP = diff(p(t),t) == DP(1,2);
disp('Costate ODEs are:');
disp([eqP1;eqP2;eqP]);
% ************************************************************************
%   Equations eqP1, eqP2 and eqP are programmed into costate_odefun
%   function. Function is implemented at end of script.
% ************************************************************************

conP1 = p1(100) == 100;
conP2 = p2(100) == 100;
conP = p(100) == 0;
disp('Costate boundary conditions are:');
disp([conP1,conP2,conP]);
% ************************************************************************
%   Boundary conditions conP1, conP2 and conP are programmed into 
%   costate_bcfun function. Implemented at end of script.
% ************************************************************************

tmesh = 0:0.1:T;

% Guess BVP solutions. Guesses are constant functions. (???) 
costate_solinit = bvpinit(tmesh,[1;1;1]);
opts = bvpset('Stats', 'on'); % enable stats for bvp solvers.

% Solve system of ODEs with 4c and 5c solvers.
disp('Solving Costate ODEs...'); disp(' ');
disp('4c solver says:');
costate_sol4c = bvp4c(@costate_odefun, @costate_bcfun, ...
    costate_solinit, opts);
% disp('5c solver says:');
% costate_sol5c = bvp5c(@costate_odefun, @costate_bcfun, ...
%   costate_solinit, opts);
costate_sol = costate_sol4c;

% Plot Costate Variables over time.
plot_solver_states(3,costate_sol,'A) Costate Variables over Time', ...
        ['p11(t)';'p22(t)';'p12(t)'],line_width);

% *** Substitute Optimal control from Ricatti Solution ***
% u_optimal = -B'*P*x = -p12(t)x1(t)-p22(t)x2(t).
% P is 2x2xM, M = length(costate_sol.x)*repelemFactor.
M = length(costate_sol.x)*repelemFactor;
P = zeros([2, 2, M]);

P(1,1,:) = repelem(costate_sol.y(1,:),repelemFactor);
P(2,2,:) = repelem(costate_sol.y(2,:),repelemFactor);
P(1,2,:) = repelem(costate_sol.y(3,:),repelemFactor);
P(2,1,:) = repelem(costate_sol.y(3,:),repelemFactor);

% Create timeseries.
dt = T/M;
x1 = zeros(1,M);
x1dot = zeros(1,M);
x2 = zeros(1,M);
x2dot = zeros(1,M);
u_optimal = zeros(1,M);

% System initial conditions
x1(1) = x10; x2(1) = x20;

for i = 1:M
    u_optimal(i) = - P(1,2,i)*x1(i) - P(2,2,i)*x2(i);
    x1dot(i) = cos(2*pi*i*dt)*x2(i);
    x2dot(i) = sin(2*pi*i*dt)*x1(i) +u_optimal(i);
    
    if i < M
        x1(i+1) = x1(i) + x1dot(i)*dt;
        x2(i+1) = x2(i) + x2dot(i)*dt;
    end
end

% Plot system results.
plot_time = (1:M)*dt;

figure;
hold on;
plot(plot_time,x1,'LineWidth',line_width);
plot(plot_time,x2,'LineWidth',line_width);
title('A) State Variables over Time');
legend('x1(t)','x2(t)');

figure;
plot(plot_time,u_optimal,'LineWidth',line_width);
title('A) Optimal Controller over Time');

%% Solution to time invariant system
disp('************************************************************');
disp('Working on time invariant problem...');
% ***State equations***
A = [0, 1; 0, 0];
AT = [0, 0; 1, 0];

% *** Ricatti Costate solutions
syms p1(t) p2(t) p(t);
P = [p1(t), p(t); p(t), p2(t)];
DP = -(P*A+AT*P- P*(B*B')*P+Q);

eqP1 = diff(p1(t),t) == DP(1,1);
eqP2 = diff(p2(t),t) == DP(2,2);
eqP = diff(p(t),t) == DP(1,2);
disp('Costate ODEs are:');
disp([eqP1;eqP2;eqP]);
% ************************************************************************
%   Equations eqP1, eqP2 and eqP are programmed into costate_odefun2
%   function. Function is implemented at end of script.
% ************************************************************************

conP1 = p1(100) == 100;
conP2 = p2(100) == 100;
conP = p(100) == 0;
disp('Costate boundary conditions are:');
disp([conP1,conP2,conP]);
% ************************************************************************
%   Boundary conditions conP1, conP2 and conP are programmed into 
%   costate_bcfun function. Implemented at end of script.
% ************************************************************************

% Guess BVP solutions. Guesses are constant functions. (???)
costate_solinit = bvpinit(tmesh,[1;1;1]);

% Solve system of ODEs with 4c and 5c solvers.
disp('Solving Costate ODEs...'); disp(' ');
disp('4c solver says:');
costate_sol4c = bvp4c(@costate_odefun2, @costate_bcfun, ...
    costate_solinit, opts);
% disp('5c solver says:');
% costate_sol5c = bvp5c(@costate_odefun2, @costate_bcfun, ...
%   costate_solinit, opts);
costate_sol = costate_sol4c;

% Plot Costate Variables over time.
plot_solver_states(3,costate_sol,'B) Costate Variables over Time', ...
        ['p11(t)';'p22(t)';'p12(t)'],line_width);

% *** Substitute Optimal control from Ricatti Solution ***
% u_optimal = -B'*P*x = -p12(t)x1(t)-p22(t)x2(t).
% P is 2x2xM, M = length(costate_sol.x)*repelemFactor.
M = length(costate_sol.x)*repelemFactor;
P = zeros([2, 2, M]);

P(1,1,:) = repelem(costate_sol.y(1,:),repelemFactor);
P(2,2,:) = repelem(costate_sol.y(2,:),repelemFactor);
P(1,2,:) = repelem(costate_sol.y(3,:),repelemFactor);
P(2,1,:) = repelem(costate_sol.y(3,:),repelemFactor);

% Create timeseries.
dt = T/M;
x1 = zeros(1,M);
x1dot = zeros(1,M);
x2 = zeros(1,M);
x2dot = zeros(1,M);
u_optimal = zeros(1,M);

% System initial conditions
x1(1) = x10; x2(1) = x20;

for i = 1:M
    u_optimal(i) = - P(1,2,i)*x1(i) - P(2,2,i)*x2(i);
    x1dot(i) = cos(2*pi*i*dt)*x2(i);
    x2dot(i) = sin(2*pi*i*dt)*x1(i) +u_optimal(i);
    
    if i < M
        x1(i+1) = x1(i) + x1dot(i)*dt;
        x2(i+1) = x2(i) + x2dot(i)*dt;
    end
end

% Plot system results.
plot_time = (1:M)*dt;

figure;
hold on;
plot(plot_time,x1,'LineWidth',line_width);
plot(plot_time,x2,'LineWidth',line_width);
title('B) State Variables over Time');
legend('x1(t)','x2(t)');

figure;
plot(plot_time,u_optimal,'LineWidth',line_width);
title('B) Optimal Controller over Time');

%% ************************************************************************
%   BVP ODE and boundary condition functions. 
% *************************************************************************
%
% A) BVP for Ricatti Costate Solution to time varying system.
%       P1 is P11, P2 is P22, P3 is P12=P21.
function dPdt = costate_odefun(t,p)
    dPdt = zeros(3,1);
    dPdt(1) = p(3)^2 - 2*sin(2*pi*t)*p(3) - 1;
    dPdt(2) = p(2)^2 - 2*cos(2*pi*t)*p(3) - 1;
    dPdt(3) = p(2)*p(3) - cos(2*pi*t)*p(1) - sin(2*pi*t)*p(2);

end

% Boundary conditions are for p(tf) so bc is function only of ptf.
function bc = costate_bcfun(p0, ptf)
    % p1(tf) = p2(tf) = 100; p12(tf) = 0;
    bc = [ptf(1) - 100; ptf(2) - 100; ptf(3)];
end

% B) BVP for Ricatti Costate Solution to time invariant system.
%       P1 is P11, P2 is P22, P3 is P12=P21.
function dPdt = costate_odefun2(t,p)
    dPdt = zeros(3,1);
    dPdt(1) = p(3)^2 - 1;
    dPdt(2) = p(2)^2 - 2*p(3) - 1;
    dPdt(3) = p(2)*p(3) - p(1);

end
%% ************************************************************************
function plot_solver_states(noStates,solution, titleStr, ...
        legendStr,line_width)
    figure();
    hold on;
    for i = 1:noStates
        plot(solution.x,solution.y(i,:),'LineWidth',line_width);
    end
    title(titleStr);
    legend(legendStr);
    xlabel('Time t')
    hold off;
end