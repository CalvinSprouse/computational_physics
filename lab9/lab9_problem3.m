% lab9_problem3

% This script simulates a rocket launched from the surface of Earth.
% Note that I am Imperial Unit System intolerant.

% author: Calvin Sprouse
% date: 2024 March 01
% class: PHYS361
% assignment: Lab 9 Problem 3

% define units for conversion
u = symunit();

% define simulation parameters
% constant thrust [N]
[F_T, ~] = separateUnits(unitConvert(8000 * u.lbf, u.N));
F_T = double(F_T);
% fuel loss rate in terms of weight [N/s]
[dwdt, ~] = separateUnits(unitConvert(80 * u.lbf, u.N));
dwdt = double(dwdt);
% initial rocket weight [N]
[w_0, ~] = separateUnits(unitConvert(3000 * u.lbf, u.N));
w_0 = double(w_0);
% Earth surface gravity [m/s^2]
g = 9.80665;
% time to solve [s]
t_span = [0, 3];

% define simulation functions
% weight of the rocket as a function of time [N]
woft = @(t) w_0 - dwdt.*t;
% drag force of the rocket as a function of velocity [N]
Dofv = @(v) 0.005.*g.*(v.^2);

% make an anonymous reference to the solver
dudt_func = @(t, u) f(t, u, F_T, g, woft, Dofv);

% solve using ode45
[t_bi, u_bi] = ode45(dudt_func, t_span, [0,0]);
x_bi = u_bi(:,1);
v_bi = u_bi(:,2);

% solve using rk4
[t_rk4, u_rk4] = rk4(dudt_func, t_span(1), t_span(2), 0.05, [0,0]);
x_rk4 = u_rk4(1,:);
v_rk4 = u_rk4(2,:);


% make comparison plots for method validation
figure
plt = plot(x_bi, v_bi, "--w");
hold on;
scatter(x_rk4, v_rk4, 24, "cyan", "filled");
hold off;
plt(1).Parent.FontWeight = "normal";
plt(1).Parent.FontName = "Times New Roman";
plt(1).Parent.FontSize = 12;
xlabel("y(t) [m]");
ylabel("v(t) [m/s]");
title("Comparing methods", FontWeight="normal");
legend("ode45", "rk4", Location="northwest", EdgeColor="none");
exportgraphics(gcf, "rocket_comp.pdf");


% make 4-panel plot for results
% find acceleration at each time step using dudt func
aoft = zeros(length(t_rk4), 1);
for i = 1:length(t_rk4)
	aatt = dudt_func(t_rk4(i), u_rk4(:,i));
	aoft(i) = aatt(2);
end
% start figure
figure
% position plot
subplot(2,2,1);
plot(t_rk4, x_rk4, "-w");
xlabel("time [s]");
ylabel("altitude [m]");
title("rocket altitude", FontWeight="normal");
% velocity plot
subplot(2,2,2);
plot(t_rk4, v_rk4, "-w");
xlabel("time [s]");
ylabel("vertical velocity [m/s]");
title("rocket velocity", FontWeight="normal");
% acceleration plot
subplot(2,2,3);
plot(t_rk4, aoft, "-w");
xlabel("time [s]");
ylabel("vertical acceleration [m/s^{2}]");
title("rocket acceleration", FontWeight="normal");
% method validation plot
subplot(2,2,4);
plot(x_bi, v_bi, "-w", x_rk4, v_rk4, ".cyan");
xlabel("altitude [m]");
ylabel("velocity [m/s]");
title('validation ode45 vs. rk4', FontWeight="normal");
exportgraphics(gcf, "rocket_panels.pdf");


% define a function to model the rocket behavior for ode45
function dudt = f(t, u, F_T, g, woft, Dofv) 
    % u(1)=x, u(2)=v

	% pre-calculation steps
	% weight as a function of time [N]
	w = woft(t);
	% drag as a function of velocity [N]
	D = Dofv(u(2));

    % set the shape of the dudt array
    dudt = zeros(2, 1);
    
    % dx/dt
    dudt(1) = u(2);
    % dv/dt
    dudt(2) = (g./w).*(F_T - w - D);
end


% define a 4th-order Runge-Kutta function to model rocket behavior
function [t, u] = rk4(f, tmin, tmax, h, u0)
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
		k1 = h.*f(t(i-1), u(:, i-1));
		k2 = h.*f(t(i-1) + (h/2), u(:, i-1) + (1./2).*k1);
		k3 = h.*f(t(i-1) + (h/2), u(:, i-1) + (1./2).*k2);
		k4 = h.*f(t(i-1) + h, u(:, i-1) + k3);

		% save values
		u(:, i) = u(:, i-1) + (1./6).*(k1 + 2.*k2 + 2.*k3 + k4);
		t(i) = t(i-1) + h;
	end
end


