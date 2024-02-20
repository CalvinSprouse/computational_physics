% lab7_problem4

% This script integrates tabular data using the trapezoid method.

% author: Calvin Sprouse
% date: 2024 February 19
% class: PHYS361
% assignment: Lab 7 Problem 4

% init workspace
clear;

% define time data [s]
t = 0:5:60;
% define raw acceleration data [m/s^2]
a = [0, 3, 8, 20, 33, 42, 40, 48, 60, 12, 8, 4, 3];
% calculate normalized acceleration data [ ]
a_norm = ( (1/9.8) .* a ).^2.5;

% calculate the HSI using TrapInt
hsi_trapint = TrapInt(t, a_norm);
% using matlab trapz
hsi_trapz = trapz(t, a_norm);
% compare
diff = abs(hsi_trapint - hsi_trapz);

% print
fprintf("HSI (my function): %g\nHSI (matlab): %g\ndifference: %g\n", hsi_trapint, hsi_trapz, diff);


% define trapezoid integrator
function area = TrapInt(x_data, y_data)
	arguments (Input)
		x_data (1,:) {mustBeNumeric, mustBeReal}
		y_data (1,:) {mustBeNumeric, mustBeReal, mustBeEqualSize(x_data, y_data)}
	end

	arguments (Output)
		area (1,1) {mustBeNumeric, mustBeReal}
	end

	% define area to sum with
	area = 0;

	% iterate over data
	for i = 1:(length(x_data)-1)
		% find the height of the current trapezoid
		% by approximating a rectangle
		h = y_data(i)/2 + y_data(i+1)/2;

		% find the width of the current rectangle
		w = x_data(i+1) - x_data(i);

		% add the new area to the area calculation
		area = area + h*w;
	end
end

% define custom validation function
function mustBeEqualSize(a,b)
    if ~isequal(size(a),size(b))
        eid = "Size:notEqual";
        msg = "Size of first input must equal size of second input.";
        error(eid, msg)
    end
end