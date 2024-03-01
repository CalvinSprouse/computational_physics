% second_order_ODE_example_tofill

% initialize workspace
clear;
close all;

% define time parameters
ti = 0;
tf = 1.2;

% call built-in solver for function
[t_func, u_func] = ode45(@f, [ti, tf], [1, 2]);

% define exact solution
% t_exact = ti:0.01:tf;
t_exact = linspace(ti, tf, 10);
x_exact = exp(-2.*t_exact) .* (8.*exp(5.*t_exact) - 3)./5;
y_exact = 2 .* exp(-2.*t_exact) .* (2 .* exp(5.*t_exact) + 3) ./ 5;

% make comparison plot
figure;
plt = plot(u_func(:,1).', u_func(:,2).', "-w", x_exact, y_exact, "xw");
plt(1).Parent.FontWeight = "normal";
plt(1).Parent.FontName = "Times New Roman";
plt(1).Parent.FontSize = 12;
xlabel("x(t)");
ylabel("y(t)");
title("Comparison of ode45 and analytic solutions to ODE", FontWeight="normal");
legend("ode45", "analytic", Location="northwest", EdgeColor="none");
exportgraphics(gcf, "ode_comp.pdf");

% write function here
function dudt = f(~, u)
	dudt = zeros(2, 1);
	dudt(1) = 2*u(1) + 2*u(2);
	dudt(2) = 2*u(1) - u(2);
end