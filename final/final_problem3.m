% final_problem3

% This script simulates a ball bounding.

% Author: Calvin Sprouse
% Date: 2024 March 14
% Assignment: PHYS 361 Final Problem 3

% init workspace
clear;

% define constants
% surface level gravitational acceleration [m/s^2]
g = 9.8;
% air density [kg/m^3]
rho = 1.2;
% drag coefficient
C_D = 0.6;
% cross sectional area [m^2]
A = 3e-3;
% mass [kg]
m = 0.06;

% calculate pre factor
pre = C_D.*A.*rho./(-2.*m);

% define initial values
% positions [m]
x_0 = 0;
y_0 = 1;
% velocities [m/s]
v_x0 = 5;
v_y0 = 0;
u_0 = [x_0, y_0, v_x0, v_y0];

% define the time span to simulate [s]
t_span = [0, 0.5];
h_step = 1e-4;

% define a reference to dudt_solver
solver = @(t, u) dudt_solver(t, u, g, pre);


% solve using built-in solver (for comprison only)
[tbi, ubi] = ode45(solver, t_span, u_0);
xbi = ubi(:, 1);
vxbi = ubi(:, 3);
ybi = ubi(:, 2);
vybi = ubi(:, 4);

% solve using rk4 solver
[trk4, urk4] = rk4(solver, t_span, h_step, u_0);
xrk4 = urk4(1, :);
vxrk4 = urk4(3, :);
yrk4 = urk4(2, :);
vyrk4 = urk4(4, :);

% find the index at which y is closest to 0
% this sets the xlimits
[y_intercept, i] = min(abs(yrk4));
x_max = xrk4(i);


% plot the comparison
% figure
% plt = plot(xrk4, yrk4, "-w", xbi, ybi, "--w");
% plt(1).Parent.FontWeight = "normal";
% plt(1).Parent.FontName = "Times New Roman";
% plt(1).Parent.FontSize = 12;
% xlabel("x [m]");
% ylabel("y [m]");
% title("Comparison of methods ode45 and rk4", FontWeight="normal");
% exportgraphics(gcf, "comparison.png");

% plot the position of the ball from rk4 until ball hits the ground
figure
plt = plot(xrk4, yrk4, "-w", LineWidth=2);
plt.Parent.FontWeight = "normal";
plt.Parent.FontName = "Times New Roman";
plt.Parent.FontSize = 12;
ylim([0, 1.1]);
xlim([0, x_max])
xlabel("x-position [m]");
ylabel("y-position [m]");
title("Trajectory of tennis ball thrown horizontally on Earth", FontWeight="normal");
exportgraphics(gcf, "rk4_trajectory.png");


% define dudt function for ode45
function dudt = dudt_solver(~, u, g, pre) 
    % u(1)=x, u(2)=y
	% u(3)=v_x, u(4)=v_y

    % set the shape of the dudt array
    dudt = zeros(4, 1);
    
    % dx/dt = v_x
    dudt(1) = u(3);
	% dy/dt = v_y
	dudt(2) = u(4);

	% calculate velocity
	v = sqrt( u(3).^2 + u(4).^2 );

    % dv_x/dt = 
    dudt(3) = pre .* v .* u(3);
	% dv_y/dt =
	dudt(4) = -g + (pre .* v .* u(4));
end

% define 4th-order runge-kutta solver using dudt solver
function [t, u] = rk4(dudt_solver, t_span, h, u0)
	% extract t values
	tmin = t_span(1);
	tmax = t_span(2);

    % define the number of steps to take
    N = round(tmax-tmin)/h;

	% define the time array
	t = zeros(1, N);

	% define u array
	u = zeros(length(u0), N);

    % define initial values
    t(1) = tmin;
    u(:, 1) = u0.';

	% solver loop
    for i=2:N
		% calculate k-terms
        k1 = h.*dudt_solver(t(i-1), u(:, i-1));
        k2 = h.*dudt_solver(t(i-1) + (h/2), u(:, i-1) + (1./2).*k1);
		k3 = h.*dudt_solver(t(i-1) + (h/2), u(:, i-1) + (1./2).*k2);
		k4 = h.*dudt_solver(t(i-1) + h, u(:, i-1) + k3);

		% save values
        u(:, i) = u(:, i-1) + (1./6).*(k1 + 2.*k2 + 2.*k3 + k4);
        t(i) = t(i-1) + h;
	end
end