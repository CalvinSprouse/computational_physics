% define total time [s]
ttot = 1.2;
% define time step [s]
h = 0.4;

% define initial conditions
x0 = 1;
y0 = 2;

% solve using matlab built in solvers
[t, u] = ode45(@f, [0, ttot], [x0; y0]);

% extract x and y values
x_solver = u(:, 1);
y_solver = u(:, 1);

% calculate exact solutions
t_exact = 0:0.01:1.2;
x_exact = (1./5) .* exp(-2.*t_exact) .* (8.*exp(5.*t_exact) - 3);
y_exact = (2./5) .* exp(-2.*t_exact) .* (2.*exp(5.*t_exact) + 3);

% plot and compare
figure;
% plt = scatter(Ioft(:,1)*1e3, Ioft(:,2)*1e6, 12, "k", "filled");
plt = plot(x_solver, y_solver, "-k");
hold on
plot(x_exact, y_exact, "--k");
% plt.Parent.Box = "on";
plt.Parent.FontWeight = "normal";
plt.Parent.FontName = "Times New Roman";
plt.Parent.FontSize = 12;
xlabel("x(t)");
ylabel("y(t)");
title("Comparison of ode45 and analytic solutions to ODE", FontWeight="normal");
legend("ode45", "analytic", Location="northwest");
exportgraphics(gcf, "ode_comp.pdf");


% define the dudt solver function
function dudt = f(~, u)
	dudt = zeros(2,1);
	dudt(1) = 2*u(1) + 2*u(2);
	dudt(2) = 2*u(1) - u(2);
end