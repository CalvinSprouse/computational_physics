% riemann_sum

% this script calculates riemann_sums.

% author: Calvin Sprouse
% date: 2024 February 11

% init workspace
clear;

% define rules
rules = ["midpoint", "left", "right", "trapezoid"];

% define function, integration limits, and rule
func = @(x) exp(x);
limits = [0, 1];
rule = rules(1);

% define the number of points [#]
N = 1e6;

% calculate the spacing Deltax and rectangle left positions
Dx = (limits(2) - limits(1)) / N;
rec_leftpos = limits(1) : Dx : limits(2)-Dx;

% calculate rectangle heights
if strcmp(rule, rules(1))
	% midpoint
	rec_heights = func(rec_leftpos + Dx/2);
elseif strcmp(rule, rules(2))
	% left
	rec_heights = func(rec_leftpos);
elseif strcmp(rule, rules(3))
	% right
	rec_heights = func(rec_leftpos + Dx);
elseif strcmp(rule, rules(4))
	% trapezoid
	rec_heights = func(rec_leftpos) + ( func(rec_leftpos + Dx) - func(rec_leftpos) ) ./ 2;
else
	error("Rule '%s' not recognized. Rule must be '%s', '%s', or '%s'.", rule, rules(1), rules(2), rules(3));
end

% calculate rectangle heights
rec_areas = rec_heights .* Dx;

% calculate sum
area = sum(rec_areas);
