% lab7_problem3

% This script implements Simpson's rule for integration.

% author: Calvin Sprouse
% date: 2024 February 19
% class: PHYS361
% assignment: Lab 7 Problem 3

% init workspace
clear;

% problem a.
% define the radius of the earth R [m]
R = 6371e3;
% define the surface level gravitational acceleration [m/s^2]
g0 = 9.81;
% define the mass of the object being raised [kg]
m = 500;

% define the gravitational acceleration function g(y) [m/s^2]
g = @(y) g0 .* (R.^2) ./ ( (R+y).^2 );

% calculate the change in potential by Simpson's rule
% from ground level to a height h [m]
h = 800e3;
DeltaU = m * SimpsonInt(g, 0, h);
DeltaU_matlab = m * integral(g, 0, h);
% print
diff = abs(DeltaU - DeltaU_matlab);
fprintf("The change in potential energy is %g Joules.\n", DeltaU);
% fprintf("The change in potential energy is %g Joules.\n", DeltaU_matlab);
fprintf("The difference between my function and MATLABs integrate is %g Joules.\n", diff);


% problem b.
% define the semi-major axis of pluto [km]
a_pluto = 5.9065e9;
% define the semi-minor axis of pluto [km]
b_pluto = 5.7208e9;
% define the period of plutos orbit [h]
T_pluto = 248 * 365 * 24;

% define a function to calculate the perimeter of an ellipse
k = sqrt( (a_pluto.^2 - b_pluto.^2) ./ a_pluto );
p = @(th) sqrt( 1 - (k.^2 .* sin(th).^2) );

% calculate the perimeter of plutos orbit by Simpson's rule
th_fin = pi/2;
P = 4*a_pluto * SimpsonInt(p, 0, th_fin);
% and by matlab
P_matlab = 4*a_pluto * integral(p, 0, th_fin);
% print
diff = abs(P - P_matlab);
fprintf("\nThe perimeter of Plutos orbit is %g kilometers.\n", P);
% fprintf("The perimeter of Plutos orbit is %g kilometers.\n", P_matlab);
fprintf("The difference between my function and MATLABs integrate is %g kilometers.\n", diff);

% calculate the average speed of Pluto
v = P / T_pluto;
v_matlab = P_matlab / T_pluto;
% print
diff = abs(v - v_matlab);
fprintf("\nThe average speed of Pluto is %g kilometers per hour.\n", v);
% fprintf("The average speed of Pluto is %g kilometers per hour.\n", v_matlab);
fprintf("The difference between my function and MATLABs integrate is %g kilometers per hour.\n", diff);


% Simpson integrator
function area = SimpsonInt(f, a, b, N)
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
