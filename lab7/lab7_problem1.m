% lab7_problem1

% This script calcualtes a numerical derivative for time-position-like data.

% author: Calvin Sprouse
% date: 2024 February 18
% class: PHYS361
% assignment: Lab 7 Problem 1

% initialize workspace
clear;

% define data
time = 0:5:40;
position = [0, 20, 53, 295, 827, 1437, 2234, 3300, 4658];

% calculate derivative
velocity = FirstDeriv(time, position);


% define functions to calculate derivatives
function deriv = FirstDeriv(time, position)
	arguments (Input)
		time (1,:) double {mustBeNumeric, mustBeReal}
		position (1,:) double {mustBeNumeric, mustBeReal, mustBeEqualSize(time, position)}
	end

	arguments (Output)
		deriv (1,:) double
	end

	% define an array of derivatives
	deriv = zeros(1, length(time));

	for i = 1:length(time)
		% check which region of points we are using
		% this determines the method
		if i <= 2
			% forward difference method
			deriv(i) = (position(i+1) - position(i)) / (time(i+1) - time(i));
		elseif length(time) - i <= 2
			% backward difference method
			deriv(i) = (position(i) - position(i-1)) / (time(i) - time(i-1));
		else
			% central difference method
			deriv(i) = (position(i+1) - position(i-1)) / (time(i+1) - time(i-1));
		end
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