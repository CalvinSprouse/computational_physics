% lab9_problem1

% This script solves a first order ODE.
% Note that colors were chosen on the dark mode screen.

% author: Calvin Sprouse
% date: 2024 March 01
% class: PHYS361
% assignment: Lab 9 Problem 1

% define inputs for ode solving
x_span = [1, 3];
y_init = 1;

% (a.) define the ODE as an anonymous function and solve with ode45
dydx_anon = @(x,y) x - x.*y./2;
[x_anon, y_anon] = ode45(dydx_anon, x_span, y_init);

% (b.) define the ODE as a function and solve with ode45
[x_func, y_func] = ode45(@f, x_span, y_init);

% (c.) define the exact solution
x_exact = 1:0.1:3;
y_exact = 2 - exp( (1 - (x_exact.^2))./4 );

% plot solutions for comparison
figure;
plt = plot(x_exact, y_exact, "-w");
hold on;
scatter(x_anon, y_anon, 28, "cyan", "filled");
scatter(x_func, y_func, 28, "magenta");
hold off;
plt(1).Parent.FontWeight = "normal";
plt(1).Parent.FontName = "Times New Roman";
plt(1).Parent.FontSize = 12;
xlabel("x");
ylabel("y");
title("Comparison of ODE solving methods", FontWeight="normal");
legend("analytic solution", "anonymous function", "defined function", Location="northwest", EdgeColor="none");
exportgraphics(gcf, "p1_ode_comp.pdf");


% define the user-defined ODE
function dydx = f(x, y)
	dydx = x - x.*y./2;
end
