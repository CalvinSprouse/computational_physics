% anemometer_analysis

% this program analyses anemometer data and fits to model
% model: windspeed = A * exp( B * voltage )
% linearized model: ln(windspeed) = ln(A) + (B * voltage)

% author: Calvin Sprouse
% date: 2024 February 01
% assignment: phys363 lab 5

% load workspace data
load("anemometer.mat");

% linearize y-component
lin_vel = log(vel);


% plot linearized data
s = scatter(volts, lin_vel, 32, "k", "filled");
% configure plot options
s.Parent.Box = "on";
s.Parent.FontWeight = "normal";
s.Parent.FontName = "Times New Roman";
s.Parent.FontSize = 16;
% specify limits
ylim([1, round(max(lin_vel))+1]);
xlim([round(min(volts)), round(max(volts))]);
% set labels
xlabel("Voltage [Volts]");
ylabel("ln(wind speed) [ln(m/s)]");
title("Linearization of anemometer data", FontWeight="normal");
% save plot
exportgraphics(gcf, "linearized_anemometer.png");


% calculate linearization coefficients through function
[slope, intercept] = linearRegression(volts, lin_vel);

% calculate A and B from wind speed model
A_fit = exp(intercept);
B_fit = slope;

% generate model data
model_vlt = round(min(volts)) : 0.01 : max(volts);
model_ws = A_fit .* exp(B_fit .* model_vlt);


% plot data in original form with model fit
figure;
s = scatter(volts, vel, 32, "k", "filled");
hold on
plot(model_vlt, model_ws, "--k");
% plot(model_vlt, model_ws_ml, "--r");
hold off
% configure plot options
s.Parent.Box = "on";
s.Parent.FontWeight = "normal";
s.Parent.FontName = "Times New Roman";
s.Parent.FontSize = 16;
% set labels
ylabel("Voltage [Volts]");
xlabel("Wind speed [m/s]");
title("Anemometer data with model fit", FontWeight="normal");
% add legend
legend("Data", "Model", Location="southeast");
% save plot
exportgraphics(gcf, "anemometer_fit.png");


% define a function to preform linear regression
function [slope, intercept] = linearRegression(x, y)
  % calculate the slope a and intercept b for some set of points
  % calculate the length of vectors
  N_x = length(x);
  N_y = length(y);

  % raise an error for N_x != N_y
  if N_x ~= N_y
	  err_size = sprintf("size(x)=%d, size(y)=%d.", N_x, N_y);
	  error("Cannot do linear regression on unequaly sized arrays: " + err_size)
  end

  % raise an error for N_x <= 1 entry
  if N_x <= 1
	  error("Cannot do linear regression with less than 2 points.")
  end

  % calculate the various sums
  S_x = sum(x);
  S_y = sum(y);
  S_xy = sum(x .* y);
  S_xx = sum(x .^ 2);

  % calculate the regression constants
  slope = ( (N_x * S_xy) - (S_x * S_y) ) / ( (N_x * S_xx) - (S_x .^ 2) );
  intercept = ( (S_xx * S_y) - (S_xy * S_x) ) / ( (N_x * S_xx) - (S_x .^ 2) );
end