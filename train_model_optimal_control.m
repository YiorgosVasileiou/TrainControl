%% ***** Optimal Train Control - Problem Parameters *****
% Given model parameters for problem simulation.
clear;
close all;

disp('******************************************************************');
global k1 k2 k3 k4 R Imin Imax  c1 c2 x1f x01 x02
k1 = 0.5; k2 = 0.1; k3 = 1; k4 = 10;
c1 = 1000; c2 = 1000; R = 0.3;
x1f = 10;
Imin = -2; Imax = 2;
T = 10;
x01 = 0; x02 = -20;

lwidth = 1;

%% ***** Dynamic system simulation w/ optimal controller. *****
% PARTS 5 and 6
tmesh = 0:0.01:T;
opts = bvpset('Stats', 'on'); % enable stats for bvp solvers.

% Guess BVP solutions for unbounded input. Guesses are constant functions.
solinit = bvpinit(tmesh,[1;1;1;1]);

% Solve system of ODEs with 4c solver.
disp('Solving State ODEs...'); disp(' ');
disp('4c solver says:');
disp('Unbound input:');
sol4c_unbound = bvp4c(@odefun_unbound, @bcfun, solinit, opts);

% Guess BVP solutions for bounded input.
% Guess is unbound input problem.
disp('Bound input:');
solinit = bvpinit(sol4c_unbound,[0 T]);
sol4c_bound = bvp4c(@odefun_bound, @bcfun, solinit, opts);

t   = sol4c_bound.x;
x1  = sol4c_bound.y(1,:);
x2  = sol4c_bound.y(2,:);
p1  = sol4c_bound.y(3,:);
p2  = sol4c_bound.y(4,:);

figure; hold on;
plot(t,x1,'LineWidth',lwidth);
plot(t,x2,'LineWidth',lwidth);
title('State Variables over Time');
legend('Position x1(t)','Speed x2(t)');

figure; hold on;
plot(t,p1,'LineWidth',lwidth);
plot(t,p2,'LineWidth',lwidth);
title('Costate Variables over Time');
legend('p1(t)','p2(t)');

% Calculated optimal controller used to drive system.
M = length(t);
u_optimal = zeros(1,M);
test = zeros(1,M);
for i = 1:M
    [u_optimal(i), ~] = opt_controller_bound(p2(i), x2(i));
end

figure;
plot(t,u_optimal,'LineWidth',lwidth);
title('Optimal Controller over Time');

disp(' ');
disp('Calculated cost for chosen control is:');
J = c1*(x1(M)-x1f)^2+c2*x2(M)^2+trapz(k4*x2.*u_optimal+R*(u_optimal.^2));
disp(J);
%% ************************************************************************
% TPBVP functions for ODEs and BCs
% *************************************************************************

% Total system s ODEs for simulation - use nonbounded input.
% s(1) = x1, s(2) = x2, s(3) = p1, s(4) = p2.
function dSdt = odefun_unbound(t,s)
    global k1 k2 k3 k4 
    [u,~] = opt_controller_unbound(s(4),s(2));
    dSdt = zeros(4,1);
    dSdt(1) = s(2);
    dSdt(2) = -k1*s(2)-k2*s(2)^2+k3*u;
    dSdt(3) = 0;
    dSdt(4) = -k4*u + k1*s(2) - s(3) + 2*k2*s(2)*s(4);
end

% Total system s ODEs for simulation - use bounded input.
% s(1) = x1, s(2) = x2, s(3) = p1, s(4) = p2.
function dSdt = odefun_bound(t,s)
    global k1 k2 k3 k4 
    [u,~] = opt_controller_bound(s(4),s(2));
    dSdt = zeros(4,1);
    dSdt(1) = s(2);
    dSdt(2) = -k1*s(2)-k2*s(2)^2+k3*u;
    dSdt(3) = 0;
    dSdt(4) = -k4*u + k1*s(2) - s(3) + 2*k2*s(2)*s(4);
end

% Boundary conditions are for x(0) and p(tf) so bc is function only of ptf.
function bc = bcfun(s0, stf)
    % x1(0) = x2(0) = 0
    % p1(tf) = 2*c1*(x1(tf)-x1f)
    % p2(tf) = 2*c2*x2(tf)
    global c1 c2 x1f x01 x02
    bc = [s0(1)-x01, ... s0(1) is x1(0)
        s0(2)-x02,  ... s0(2) is x2(0)
        stf(3)-2*c1*(stf(1)-x1f),... stf(3) is p1(tf), stf(1) is x1(tf)  
        stf(4)-2*c2*stf(2)];    % stf(4) is p2(tf), stf(2) is p2(tf)
end

%% ***** Optimal Controller Functions *****
function [u,t] = opt_controller_bound(p2, x2)
    global k3 k4 R Imax Imin
    t = k3*p2 + k4*x2;
    
    u = -t/(2*R);
    if u < Imin
        u = Imin;
    elseif u > Imax
         u = Imax;
    end
end

function [u,t] = opt_controller_unbound(p2, x2)
    global k3 k4 R
    t = k3*p2 + k4*x2;
    
    u = -t/(2*R);
end

%% ***** Linearizarion around optimal curve *****
% PART 7

% Linear model for error near optimal curve.
% y1dot = y2;
% y2dot = -k1 - 2k2x2 + k3u
% Costate dynamic equations - Ricatti solutions.
% 
