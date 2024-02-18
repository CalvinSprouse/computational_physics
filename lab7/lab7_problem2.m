% lab7_problem2

% This script calcualtes two numerical derivatives for time-position-like data.

% author: Calvin Sprouse
% date: 2024 February 18
% class: PHYS361
% assignment: Lab 7 Problem 2

% initialize workspace
clear;

% define data
time = 0:10:120;
position = [-8, 241, 1244, 2872, 5377, 8130, 11617, 15380, 19872, 25608, 31412, 38309, 44726];

% calculate derivatives
velocity = FirstDeriv(time, position);
acceleration = FirstDeriv(time, velocity);

% make plot of position vs time
plt = scatter(time, position, 32, "filled", "ok");
plt.Parent.Box = "on";
plt.Parent.FontWeight = "normal";
plt.Parent.FontName = "Times New Roman";
plt.Parent.FontSize = 16;
xlabel("Time [s]");
ylabel("Position [m]");
title("Position of an object", FontWeight="normal");
exportgraphics(gcf, "position.pdf");

% make plot of velocity vs time
plt = scatter(time, velocity, 32, "filled", "ok");
plt.Parent.Box = "on";
plt.Parent.FontWeight = "normal";
plt.Parent.FontName = "Times New Roman";
plt.Parent.FontSize = 16;
xlabel("Time [s]");
ylabel("Velocity [m/s]");
title("Velocity of an object", FontWeight="normal");
exportgraphics(gcf, "velocity.pdf");

% make plot of acceleration vs time
plt = scatter(time, acceleration, 32, "filled", "ok");
plt.Parent.Box = "on";
plt.Parent.FontWeight = "normal";
plt.Parent.FontName = "Times New Roman";
plt.Parent.FontSize = 16;
xlabel("Time [s]");
ylabel("Acceleration [m/s^2]");
title("Acceleration of an object", FontWeight="normal");
exportgraphics(gcf, "acceleration.pdf");


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