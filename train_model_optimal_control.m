%% ***** Optimal Train Control - Problem Parameters *****
% Given model parameters for problem simulation.
clear;
%close all;

disp('******************************************************************');
global k1 k2 k3 k4 R Imin Imax  c1 c2 x1f x01 x02
global x01_ic_error x02_ic_error
k1 = 0.5; k2 = 0.1; k3 = 1; k4 = 10;
c1 = 1000; c2 = 1000; R = 0.3;
x1f = 10;
Imin = -2; Imax = 2;
T = 10;
x01 = 0; x02 = 0;
x01_ic_error = 0.1; x02_ic_error = 0.3;

lwidth = 1;

%% ***** Dynamic system simulation w/ optimal controller. *****
% PARTS 5 - Find optimal control for minimizing given cost J.
% Control system dynamic equations:
% x1dot = x2
% x2dot = -k1*x2 - k2*x2^2 + k3* u_opt
% Costate dynamic equations.
% p1dot = 0
% p2dot = -k4*u + k1*x2 - p1 + 2*k2*x2*p2

tmesh = 0:0.01:T;
opts = bvpset('Stats', 'on'); % enable stats for bvp solvers.

% Guess BVP solutions for unbounded input. Guesses are constants.
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

for i = 1:M
    [u_optimal(i), ~] = opt_controller_bound(p2(i), x2(i));
end

figure;
plot(t,u_optimal,'LineWidth',lwidth);
title('Optimal Controller over Time');

disp(' ');
disp('Calculated cost for optimal control is:');
J=c1*(x1(M)-x1f)^2+c2*x2(M)^2+trapz(k4*x2.*u_optimal+R*(u_optimal.^2));
disp(J);

% **********************************************************************
% PART 6 - Drive train with previous control from different starting
% positions.

x1_err = zeros(1,M); x1dot_err = zeros(1,M);
x2_err = zeros(1,M); x2dot_err = zeros(1,M);

x1_err(1) = x01 + x01_ic_error; x2_err(1) = x02 + x02_ic_error;
dt = T/M;

% Simulate system response.
for i = 1:M
    x1dot_err(i) = x2_err(i);
    x2dot_err(i) = -k1*x2_err(i)-k2*x2_err(i)^2 + k3*u_optimal(i);
    if i < M
        x1_err(i+1) = x1_err(i) + x1dot_err(i)*dt;
        x2_err(i+1) = x2_err(i) + x2dot_err(i)*dt;
    end
end

figure; hold on;
plot(t,x1_err,'LineWidth',lwidth);
plot(t,x2_err,'LineWidth',lwidth);
title('State Variables over Time, Error in IC');
legend('Position x1(t)','Speed x2(t)');

disp(' ');
disp('Calculated cost for sub_optimal control is:');
J_err = c1*(x1_err(M)-x1f)^2+c2*x2_err(M)^2 + ...
    trapz(k4*x2_err.*u_optimal+R*(u_optimal.^2));
disp(J_err);


%% ***** Linearizarion around optimal curve *****
% PART 7 - Minimize cost J2 for linearization around optimal trajectory.

% Linear model for error near optimal curve -> LINEAR PROBLEM:
% y1dot = y2;
% y2dot = (-k1 - 2*k2*x2)y2 + k3*v
% Costate dynamic equations -> RICATTI PROBLEM:
% a = -k1 - 2*k2*x2
% p1dot = k3^2*p^2 - 2
% p2dot = k3^2*p2^2 - 2 - 2*(p+a*p2)
% ldot  = k3^2*p2*l - a*p - p1

% Try linearization around polynomial estimation of x2.
global poly_coeff
poly_coeff = polyfit(t,x2,8); 
disp(' ');
disp('Linearizing around optimal trajectory...'); disp(' ');
disp('4c solver says:');
disp('Linear problem solution:')

% Solve linearization system for error correction.
tmesh = linspace(0,T,M);
opts = bvpset('NMax',M,'Stats', 'on','Abstol',1e-05);
solinit = bvpinit(tmesh,@guess);
sol4c_linprob_ricatti = bvp5c(@odefun_linprob_ricatti, ...
    @bcfun_linprob_ricatti, solinit, opts);

t = sol4c_linprob_ricatti.x;
p11 = sol4c_linprob_ricatti.y(1,:);
p22 = sol4c_linprob_ricatti.y(2,:);
p12 = sol4c_linprob_ricatti.y(3,:);

% Calculated optimal controller used to drive linear error system.
v_optimal = zeros(1,M);
y1 = zeros(1,M); y1dot = zeros(1,M);
y2 = zeros(1,M); y2dot = zeros(1,M);

y1(1) = x01_ic_error - x01; y2(1) = x02_ic_error - x02 ;
dt = T/M;

% Simulate system response.
for i = 1:M
    v_optimal(i) = v_optimal_controller(y1(i),y2(i),p12(i),p22(i));
    a = -k1 - 2*k2*x2(i);
    
    y1dot(i) = y2(i);
    y2dot(i) = a*y2(i) + k3*v_optimal(i);
    if i < M
        y1(i+1) = y1(i) + y1dot(i)*dt;
        y2(i+1) = y2(i) + y2dot(i)*dt;
    end
end

figure; hold on;
plot(t,y1,'LineWidth',lwidth);
plot(t,y2,'LineWidth',lwidth);
title('Lin. State Variables over Time');
legend('Position Error y1(t)','Speed Error y2(t)');

figure;
plot(t,v_optimal,'LineWidth',lwidth);
title('Optimal Lin. Controller over Time');

% **********************************************************************
% PART 8 - Drive train with corrected control from different starting
% positions.

u_modified = u_optimal + v_optimal;

x1_err_mod = zeros(1,M); x1dot_err_mod = zeros(1,M);
x2_err_mod = zeros(1,M); x2dot_err_mod = zeros(1,M);

x1_err_mod(1) = x1_err(1); x2_err_mod(1) = x2_err(1);
dt = T/M;

for i = 1:M
    x1dot_err_mod(i) = x2_err_mod(i);
    x2dot_err_mod(i) = -k1*x2_err_mod(i)-k2*x2_err_mod(i)^2 + ...
        k3*u_modified(i);
    if i < M
        x1_err_mod(i+1) = x1_err_mod(i) + x1dot_err_mod(i)*dt;
        x2_err_mod(i+1) = x2_err_mod(i) + x2dot_err_mod(i)*dt;
    end
end

figure; hold on;
plot(t,u_optimal,'b--','LineWidth',lwidth);
plot(t,v_optimal,'r--','LineWidth',lwidth);
plot(t,u_modified,'m','LineWidth',lwidth);
title('Mod. Controller over Time');
legend('Optimal Controller','Error Correction', 'Modified Controller');

figure; hold on;
plot(t,x1_err_mod,'LineWidth',lwidth);
plot(t,x2_err_mod,'LineWidth',lwidth);
title('Mod. State Variables over Time, Error in IC');
legend('Position x1(t)','Speed x2(t)');

disp(' ');
disp('Calculated cost for corrected control is:');
J_err_mod = c1*(x1_err_mod(M)-x1f)^2+c2*x2_err_mod(M)^2 + ...
    trapz(k4*x2_err_mod.*u_modified+R*(u_modified.^2));
disp(J_err_mod);

%% **********************************************************************
% TPBVP functions for ODEs and BCs
% ***********************************************************************

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

% Boundary conditions are for x(0) and p(tf).
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

%% ***** Non-linear System Optimal Controller Functions *****
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

%% **********************************************************************
% TPBVP functions for ODEs and BCs - LINEAR PROBLEM
% ***********************************************************************

% Linear ODEs for error. Ricatti problem
% y(1) = p11, y(2) = p22, y(3) = p12 = p21.
function dYdt = odefun_linprob_ricatti(t,y)
    global k1 k2 k3 poly_coeff
    x2 = polyval(poly_coeff,t);
    a = - (k1 + 2*k2*x2);
    K = k3^2;
    
    dYdt = zeros(3,1);
    % Ricatti variables.
    dYdt(1) = K*y(3)^2 - 2;
    dYdt(2) = -2-2*(y(3)+a*y(2))+K*y(2)^2;
    dYdt(3) = K*y(2)*y(3)- a*y(3) - y(1);
end

% Boundary conditions for y(0) and y(tf).
function bc = bcfun_linprob_ricatti(y0, ytf)
    % p11(tf) = p22(tf) = 20
    % p12(tf) = 0
    bc = [ytf(1)-20, ytf(2)-20, ytf(3)];
end

%% ***** Linearized System Error Optimal Controller Functions *****
function v = v_optimal_controller(y1, y2, p2, p12)
    global k3
    v = -k3*(p12*y1+p2*y2);
end

function g = guess(x)
   g = [x;1;1];
end