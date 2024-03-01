% SHO_problem_to_fill

% This is a template function to write a simple harmonic oscillator solver.

% initialize workspace
clear;
close all;

% initialize constants
% define time parameters [s]
tmin = 0;
tmax = 500;
% define omega [rad/s]
omega = 0.1;

% define initial values
% initial position [m]
x0 = 0;
% initial velocity [m/s]
v0 = 1;

% package the time span and the intial values as arrays
tspan = [tmin, tmax];
initvals = [x0, v0];

% solve using built-in solver
[tbi, ubi] = ode45(@(t,u) sho(t,u,omega), tspan, initvals);
xbi = ubi(:, 1);
vbi = ubi(:, 2);

% solve using 2nd order Runge-Kutta
[trk2, urk2] = rk2(@(t,u) sho(t,u,omega), tmin, tmax, 0.1, initvals);
xrk2 = urk2(1, :);
vrk2 = urk2(2, :);

% solve using 4th order Runge-Kutta
[trk4, urk4] = rk4(@(t,u) sho(t,u,omega), tmin, tmax, 0.1, initvals);
xrk4 = urk4(1, :);
vrk4 = urk4(2, :);

% plot the comparison
figure;
plt = plot(xbi, vbi, "-w", xrk2, vrk2, "--r", xrk4, vrk4, "-b");
% plt = plot(xbi, vbi, "-w");
% plt = plot(xrk2, vrk2, "--w");
plt(1).Parent.FontWeight = "normal";
plt(1).Parent.FontName = "Times New Roman";
plt(1).Parent.FontSize = 12;
xlabel("x(t) [m]");
ylabel("v(t) [m/s]");
title("Comparing methods", FontWeight="normal");
legend("ode45", "rk2", "rk4", Location="northwest", EdgeColor="none");
exportgraphics(gcf, "ode_comp.pdf");


% make the subplots
% h=1
[t1rk4, u1rk4] = rk4(@(t,u) sho(t,u,omega), tmin, tmax, 1, initvals);
[t1rk2, u1rk2] = rk2(@(t,u) sho(t,u,omega), tmin, tmax, 1, initvals);
% h=5
[t2rk4, u2rk4] = rk4(@(t,u) sho(t,u,omega), tmin, tmax, 5, initvals);
[t2rk2, u2rk2] = rk2(@(t,u) sho(t,u,omega), tmin, tmax, 5, initvals);
% h=10
[t3rk4, u3rk4] = rk4(@(t,u) sho(t,u,omega), tmin, tmax, 10, initvals);
[t3rk2, u3rk2] = rk2(@(t,u) sho(t,u,omega), tmin, tmax, 10, initvals);
% h=50
[t4rk4, u4rk4] = rk4(@(t,u) sho(t,u,omega), tmin, tmax, 50, initvals);
[t4rk2, u4rk2] = rk2(@(t,u) sho(t,u,omega), tmin, tmax, 50, initvals);

% extract the x positions and velocities
% from rk2
x1rk2 = u1rk2(1,:);
v1rk2 = u1rk2(2,:);
x2rk2 = u2rk2(1,:);
v2rk2 = u2rk2(2,:);
x3rk2 = u3rk2(1,:);
v3rk2 = u3rk2(2,:);
x4rk2 = u4rk2(1,:);
v4rk2 = u4rk2(2,:);
% from rk4
x1rk4 = u1rk4(1,:);
v1rk4 = u1rk4(2,:);
x2rk4 = u2rk4(1,:);
v2rk4 = u2rk4(2,:);
x3rk4 = u3rk4(1,:);
v3rk4 = u3rk4(2,:);
x4rk4 = u4rk4(1,:);
v4rk4 = u4rk4(2,:);

% plot results
subplot(2,2,1);
plot(x1rk4, v1rk4, "-r", x1rk2, v1rk2, "-b", xbi, vbi, "-w");
% plot(x1rk4, v1rk4, "-r", x1rk2, v1rk2, "-b");
xlabel("x(t) [m]");
ylabel("v(t) [m/s]");
title('h=1');

subplot(2,2,2);
plot(x2rk4, v2rk4, "-r", x2rk2, v2rk2, "-b", xbi, vbi, "-w");
% plot(x2rk4, v2rk4, "-r", x2rk2, v2rk2, "-b");
xlabel("x(t) [m]");
ylabel("v(t) [m/s]");
title('h=5');

subplot(2,2,3);
plot(x3rk4, v3rk4, "-r", x3rk2, v3rk2, "-b", xbi, vbi, "-w");
% plot(x3rk4, v3rk4, "-r", x3rk2, v3rk2, "-b");
xlabel("x(t) [m]");
ylabel("v(t) [m/s]");
title('h=10');

subplot(2,2,4);
plot(x4rk4, v4rk4, "-r", x4rk2, v4rk2, "-b", xbi, vbi, "-w");
% plot(x4rk4, v4rk4, "-r", x4rk2, v4rk2, "-b");
xlabel("x(t) [m]");
ylabel("v(t) [m/s]");
title('h=50');
exportgraphics(gcf, "ode_corners.pdf");


% define dudt function for ode45
function dudt = sho(~, u, omega) 
    % u(1) is x
    % u(2) is v

    % set the shape of the dudt array
    dudt = zeros(2, 1);
    
    % dx/dt
    dudt(1) = u(2);
    % dv/dt
    dudt(2) = -(omega.^2) * u(1); 
end


% define 2nd-order Runge-Kutta solver
function [t, u] = rk2(sho, tmin, tmax, h, u0)
    % define the number of steps to take
    N = round(tmax-tmin)/h;

	% define the time array
	t = zeros(1, N);

	% define u array
	u = zeros(2, N);

    % define initial values
    t(1) = tmin;
    u(:, 1) = u0.';

	% solver loop
    for i=2:N
		% calculate k-terms
        k1 = (h./2).*sho(t(i-1), u(:, i-1));
        k2 = h.*sho(t(i-1) + (h/2), u(:, i-1) + k1);
	
		% save values
        u(:, i) = u(:, i-1) + k2;
        t(i) = t(i-1) + h;
	end
end

% define 4th-order runge-kutta
function [t, u] = rk4(sho, tmin, tmax, h, u0)
    % define the number of steps to take
    N = round(tmax-tmin)/h;

	% define the time array
	t = zeros(1, N);

	% define u array
	u = zeros(2, N);

    % define initial values
    t(1) = tmin;
    u(:, 1) = u0.';

	% solver loop
    for i=2:N
		% calculate k-terms
        k1 = h.*sho(t(i-1), u(:, i-1));
        k2 = h.*sho(t(i-1) + (h/2), u(:, i-1) + (1./2).*k1);
		k3 = h.*sho(t(i-1) + (h/2), u(:, i-1) + (1./2).*k2);
		k4 = h.*sho(t(i-1) + h, u(:, i-1) + k3);

		% save values
        u(:, i) = u(:, i-1) + (1./6).*(k1 + 2.*k2 + 2.*k3 + k4);
        t(i) = t(i-1) + h;
	end
end
