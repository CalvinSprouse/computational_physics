% lab8_problem2

% This script uses 2nd-order Runge-Kutta to model a circuit.

% author: Calvin Sprouse
% date: 2024 February 25
% class: PHYS361 Lab 8 Problem 2

% init workspace
clear;

% define circuit component values
% DC voltage [V]
V0 = 500;
% inductor [Henries]
L = 15;
% resistance function [Amps] -> [Ohms]
R = @(I) 500 + 250.*(I.^2);

% define total time and time step
% total time [s]
ttot = 0.1;
% time step [s]
Deltat = 0.005;

% define the dI/dt function and initial condition
% [Amps/s]
dIdt = @(t,I) (V0/L) - (1/L)*R(I)*I;
% [Amps]
I0 = 0;

% apply the runge_kutta function to dIdt
Ioft = runge_kutta_2(dIdt, 0, ttot, Deltat, I0);

% plot I(t)
figure;
plt = scatter(Ioft(:,1), Ioft(:,2), 36, "k", "filled");
plt.Parent.Box = "on";
plt.Parent.FontWeight = "normal";
plt.Parent.FontName = "Times New Roman";
plt.Parent.FontSize = 16;
xlabel("time [s]");
ylabel("Current [Amps]");
title("Circuit current by 2nd-order Runge-Kutta", FontWeight="normal");
exportgraphics(gcf, "current_time.pdf");


% define a second order runge-kutta solver function
function func_vals = runge_kutta_2(func, a, b, Delta, fa)
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

		% calculate k1, k2
		k1 = func(t, func_now);
		k2 = func(t + Delta/2, func_now + (Delta/2)*k1);

		% calculate next func val
		func_next = func_now + Delta*k2;

		% save to array
		func_vals(i+1,:) = [t+Delta, func_next];
	end
end