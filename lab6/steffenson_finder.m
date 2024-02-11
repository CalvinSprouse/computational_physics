% steffenson_func

% author: Calvin Sprouse
% date: 2024 February 10
% problem: 3

function [x, err] = steffenson_finder(func, xi, iter_max, err_thresh)
	arguments
		func (1,1) function_handle
		xi (1,1) double {mustBeReal}
		iter_max (1,1) double {mustBePositive,mustBeReal} = 100
		err_thresh (1,1) double {mustBePositive,mustBeReal} = 1e-6
	end

	% given a function and initial guess this applies the steffenson root
	% finding method and returns the guess (x) and approximate error (err)

	% calculate initial iteration parameters
	xip = xi;
	err = err_thresh + 1;

	% enter while loop
	iter = 0;
	while (err > err_thresh)
		% calculate next guess and error
		xip = xi - ( (func(xi).^2) / (func(xi + func(xi)) - func(xi)) );
		err = abs( (xip - xi) / xip );

		% iterate and check iteration condition
		xi = xip;
		iter = iter + 1;
		if iter > iter_max
			% raise error
			error("Iteration maximum (%d) reached.", iter_max);
		end
	end

	% return values by assignment
	x = xip;
end