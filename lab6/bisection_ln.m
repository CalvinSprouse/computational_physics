% bisection_ln

% this script calculates the natural log of a number p by solving
% the equation e^x - p = 0 using the bisection method.

% author: Calvin Sprouse
% date: 2024 February 10
% problem: 1

% initialize workspace
clear;

% define the number to find the natural log of
% ln(510) = 6.2344
p = 510;

% check input for validity
if p <= 0
	% raise error
	error("Cannot find natural log for p <= 0.");
end

% define the initial guessing values
err_approx = 1;
xr_prev = 0;

% define the threshold accuracy
err_thresh = 1e-4;

% define iteration parameters
iter = 0;
iter_max = 1e4;

% define the initial guess bounds
if p > exp(1)
	xl = exp(0);
	xu = p;
elseif p < exp(1)
	xl = -1/p;
	xu = exp(0);
else
	xr = 1;
	err_approx = 0;
end

% iterate with bisection method
while (err_approx > err_thresh) && (iter < iter_max)
	% calculate the mid point
	xr = (xu + xl) / 2;

	% determine which side to move
	if (exp(xr)-p)*(exp(xu)-p) < 0
		% replace lower bound
		xl = xr;
	else
		% replace upper bound
		xu = xr;
	end

	% calculate approximate error
	err_approx = abs( (xr - xr_prev) / xr );

	% store previous guess
	xr_prev = xr;

	% iterate
	iter = iter + 1;
end

% calculate absolute error
err_abs = abs( (log(p) - xr) / log(p) );

% output results
output = "ln(%g) ~= %g.\nSolution found after %d iterations.\nApproximate error = %g%%, absolute error = %g%%.\n";
fprintf(output, p, xr, iter, err_approx*100, err_abs*100)
