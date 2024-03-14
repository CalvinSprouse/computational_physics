% final_problem2

% This script models the heat capacity of a solid.

% Author: Calvin Sprouse
% Date: 2024 March 14
% Assignment: PHYS 361 Final Exam Problem 2

% init workspace
clear;

% define constants
% volume [cm^3]
V = 1000;
% number density [m^-3]
rho = 6.022e28;
% Debye temperature [K]
theta_D = 428;
% Boltzmann's constant [J/K]
k_B = 1.380649e-23;

% define the integrand and bounds
f = @(x) ( (x.^4).*exp(x) ./ ((exp(x) - 1).^2) );
a = 1e-6;
b = @(T) theta_D ./ T;

% define T values for plot [K]
T_vals = 5:5:500;

% calculate the heat capacity at each temperature
C_V_vals = zeros(length(T_vals), 1);
for i = 1:length(T_vals)
	% starting with a pre-factor
	C_V_pre = 9.*V.*rho.*k_B.*(T_vals(i)./theta_D).^3;

	% calculate C_V value
	C_V_vals(i) = simpson_int(f, a, b(T_vals(i)));
end

% plot C_V vs T
plt = scatter(T_vals, C_V_vals, 24, "k", "filled");
plt.Parent.Box = "on";
plt.Parent.FontWeight = "normal";
plt.Parent.FontName = "Times New Roman";
plt.Parent.FontSize = 12;
xlabel("Temperature [Kelvin]");
ylabel("Heat Capacity");
title("Heat capacity dependence on temperature for a solid", FontWeight="normal");
exportgraphics(gcf, "heat_capacity.png");


% define the Simpson integrator function
function area = simpson_int(f, a, b, N)
	arguments (Input)
		f function_handle
		a (1,1) double {mustBeNumeric, mustBeReal}
		b (1,1) double {mustBeNumeric, mustBeReal}
		N (1,1) double {mustBeNumeric, mustBeReal} = 1000
	end

	arguments (Output)
		area (1,1) double {mustBeNumeric}
	end

	% define h
	h = (b-a)/N;

	% create a vector of i values
	i = 0:2:(N-2);

	% calculate individual points to sum
	n_points = f(a + i.*h) + 4.*f(a + (i+1).*h) + f(a + (i+2).*h);

	% calculate area
	area = (h/3) * sum(n_points);
end