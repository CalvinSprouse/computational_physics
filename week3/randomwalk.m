% define distance goal [m]
dist_goal = 100;

% define an array to store particle position
% [x, y, dist]
pos = zeros(1, 3);

% define a step counter
step = 1;

% begin simulation
while pos(end, 3) < dist_goal
	% generate random direction
	theta = rand(1) * 2 * pi;

	% calculate x and y displacement
	x_disp = cos(theta);
	y_disp = sin(theta);

	% add position to array
	pos(step+1, 1) = pos(step, 1) + x_disp;
	pos(step+1, 2) = pos(step, 2) + y_disp;

	% calculate distance and append to array
	dist = sqrt(pos(step+1, 1).^2 + pos(step+1, 2).^2);
	pos(step+1, 3) = dist;

	% iterate step
	step = step + 1;
end