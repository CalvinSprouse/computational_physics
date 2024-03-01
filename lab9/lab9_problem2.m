% lab9_problem2

% This script solves a second order ODE.
% Note that colors were chosen on the dark mode screen.

%{
tic-toc built-in solver: 0.0533
tic-toc rk2: 0.0110
tic-toc rk4: 0.0188

I anticipated the built in solver would be fastest. Its possible, even
likely, that the built in solver is simply taking more steps than the rk2
and rk4 methods. Also likely that the "extra magic" beyond just a 4th-order
Runge-Kutta adds to computation time. Unsurprisingly rk2 is faster than rk4
corresponding to the fact that rk4 has just slightly more calculations.
Another possible time-save is that the rk2 and rk4 methods are not calling
a function repeatedly, which is computationally expensive due to
transfering variables, but instead remain in their function.
%}

% author: Calvin Sprouse
% date: 2024 March 01
% class: PHYS361
% assignment: Lab 9 Problem 2

% define constants
% time range [s]
t_span = [0, 500];
% angular velocity [rad/s]
omega = 0.1;
% dampening coefficient [1/s]
gamma = 0.02;

% define initial values
% initial position [m]
x0 = 0;
% initial velocity [m/s]
v0 = 1;
initvals = [x0, v0];

% solve using built-in solver
tic;
[tbi, ubi] = ode45(@(t,u) sho(t,u,omega), t_span, initvals);
toc_bi = toc;
xbi = ubi(:, 1);
vbi = ubi(:, 2);

% solve using 2nd order Runge-Kutta
tic;
[trk2, urk2] = rk2(@(t,u) sho(t,u,omega), t_span(1), t_span(2), 0.1, initvals);
toc_rk2 = toc;
xrk2 = urk2(1, :);
vrk2 = urk2(2, :);

% solve using 4th order Runge-Kutta
tic;
[trk4, urk4] = rk4(@(t,u) sho(t,u,omega), t_span(1), t_span(2), 0.1, initvals);
toc_rk4 = toc;
xrk4 = urk4(1, :);
vrk4 = urk4(2, :);

% solve the damped system using ode45
tic;
[t_damp, u_damp] = ode45(@(t,u) shodamped(t,u,omega,gamma), t_span, initvals);
toc_damp = toc;
x_damp = u_damp(:, 1);
v_damp = u_damp(:, 2);


% plot the comparison of undamped functions for rk2 and rk4 method validation
figure
plt = plot(xbi, vbi, "-w", xrk2, vrk2, "--r", xrk4, vrk4, "-b");
plt(1).Parent.FontWeight = "normal";
plt(1).Parent.FontName = "Times New Roman";
plt(1).Parent.FontSize = 12;
xlabel("x(t) [m]");
ylabel("v(t) [m/s]");
title("Comparing methods", FontWeight="normal");
legend("ode45", "rk2", "rk4", Location="northwest", EdgeColor="none");
exportgraphics(gcf, "ode_comp.pdf");


% plot the comparison of undamped and damped functions
figure
plt = plot(xbi, vbi, "--white", x_damp, v_damp, "-cyan");
plt(1).Parent.FontWeight = "normal";
plt(1).Parent.FontName = "Times New Roman";
plt(1).Parent.FontSize = 12;
xlabel("x(t) [m]");
ylabel("v(t) [m/s]");
title("Comparison of damped and undamped oscillator", FontWeight="normal");
legend("undamped", "damped", Location="northwest", EdgeColor="none");
exportgraphics(gcf, "damp_undamp.pdf");


% make the comparison subplots
tmin = t_span(1);
tmax = t_span(2);
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
x1rk2 = u1rk2(1,:); v1rk2 = u1rk2(2,:);
x2rk2 = u2rk2(1,:); v2rk2 = u2rk2(2,:);
x3rk2 = u3rk2(1,:); v3rk2 = u3rk2(2,:);
x4rk2 = u4rk2(1,:); v4rk2 = u4rk2(2,:);
% from rk4
x1rk4 = u1rk4(1,:); v1rk4 = u1rk4(2,:);
x2rk4 = u2rk4(1,:); v2rk4 = u2rk4(2,:);
x3rk4 = u3rk4(1,:); v3rk4 = u3rk4(2,:);
x4rk4 = u4rk4(1,:); v4rk4 = u4rk4(2,:);

% plot results
figure
subplot(2,2,1);
plot(x1rk4, v1rk4, "-r", x1rk2, v1rk2, "-b", xbi, vbi, "-w");
xlabel("x(t) [m]");
ylabel("v(t) [m/s]");
title('h=1');

subplot(2,2,2);
plot(x2rk4, v2rk4, "-r", x2rk2, v2rk2, "-b", xbi, vbi, "-w");
xlabel("x(t) [m]");
ylabel("v(t) [m/s]");
title('h=5');

subplot(2,2,3);
plot(x3rk4, v3rk4, "-r", x3rk2, v3rk2, "-b", xbi, vbi, "-w");
xlabel("x(t) [m]");
ylabel("v(t) [m/s]");
title('h=10');

subplot(2,2,4);
plot(x4rk4, v4rk4, "-r", x4rk2, v4rk2, "-b", xbi, vbi, "-w");
xlabel("x(t) [m]");
ylabel("v(t) [m/s]");
title('h=50');
exportgraphics(gcf, "ode_corners.pdf");


% define dudt function for ode45
function dudt = sho(~, u, omega) 
    % u(1)=x, u(2)=v

    % set the shape of the dudt array
    dudt = zeros(2, 1);
    
    % dx/dt
    dudt(1) = u(2);
    % dv/dt
    dudt(2) = -(omega.^2) * u(1); 
end


% define damped oscillator
function dudt = shodamped(~, u, omega, gamma)
	% u(1)=x, u(2)=v
	dudt = zeros(2, 1);
	% dx/dt
	dudt(1) = u(2);
	% dv/dt
	dudt(2) = -gamma*u(2) - (omega.^2)*u(1);
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

% define 4th-order runge-kutta solver
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
