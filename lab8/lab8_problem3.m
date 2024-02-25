% lab8_problem3

% This script uses a 4th-order Runge-Kutta algorithm to model current.

% author: Calvin Sprouse
% date: 2024 Februrary 25
% class: PHYS361 Lab 8 Problem 3

% init workspace
clear;

% define circuit component values
% inductor [Henries]
L = 15;
% resistance [Ohms]
R = 1000;
% frequency [Hertz]
f = 100e3;

% define the dI/dt function and initial condition
% [Amps/s]
dIdt = @(t,I) (10./L).*sin(2.*pi.*f.*t) - (R./L).*I;
% [Amps]
I0 = 0;

% define total time and time step
% total time [s]
ttot = 1e-4;
% time step [s]
Deltat = 1e-9;

% apply the runge_kutta function to dIdt
Ioft = runge_kutta_4(dIdt, 0, ttot, Deltat, I0);

% plot I(t)
figure;
plt = scatter(Ioft(:,1)*1e3, Ioft(:,2)*1e6, 12, "k", "filled");
plt.Parent.Box = "on";
plt.Parent.FontWeight = "normal";
plt.Parent.FontName = "Times New Roman";
plt.Parent.FontSize = 16;
xlabel("time [\mus]");
ylabel("Current [\muAmps]");
title("Circuit current by 4th-order Runge-Kutta", FontWeight="normal");
exportgraphics(gcf, "current_time.pdf");


% define the 4th-order Runge-Kutta solver
function func_vals = runge_kutta_4(func, a, b, Delta, fa)
	arguments (Input)
		func (1,1) function_handle
		a (1,1) {mustBeReal, mustBeNumeric}
		b (1,1) {mustBeReal, mustBeNumeric}
		Delta (1,1) {mustBeReal, mustBeNumeric}
		fa (1,1) {mustBeReal, mustBeNumeric}
	end

	arguments (Output)
		func_vals (:,2) {mustBeReal, mustBeNumeric}
	end

	% calculate the number of steps
	steps = (b - a) / Delta;

	% ensure steps are an integer value
	if round(steps) ~= steps
		error("Non-integer number of steps.")
	end

	% track values at each step
	func_vals = zeros(steps+1, 2);
	func_vals(1,:) = [a, fa];


	% iterate over steps
	for i = 1:steps
		% calculate the current input value
		t = a + (i-1)*Delta;

		% retrieve the current function value
		func_now = func_vals(i, 2);

		% calculate k1, k2, k3, k4
		k1 = func(t, func_now);
		k2 = func(t + Delta/2, func_now + (Delta/2)*k1);
		k3 = func(t + Delta/2, func_now + (Delta/2)*k2);
		k4 = func(t + Delta, func_now + Delta*k3);

		% calculate next func val
		func_next = func_now + (Delta/6)*(k1 + 2*k2 + 2*k3 + k4);

		% save to array
		func_vals(i+1,:) = [t+Delta, func_next];
	end
end