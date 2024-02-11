% euler radioactive decay

% this script solves a radioactive decay ode using Euler's method.

% author: Calvin Sprouse
% date: 2024 February 11

% init workspace
clear;

% define gamma and ttot for each sim
gamma = 0.01;
ttot = 500;

% run euler_decay_ode at Dt={1,5,10,20}
[N1, t1] = euler_decay_ode(gamma, ttot, 1);
[N5, t5] = euler_decay_ode(gamma, ttot, 5);
[N10, t10] = euler_decay_ode(gamma, ttot, 10); 
[N20, t20] = euler_decay_ode(gamma, ttot, 20);

% plot results
subplot(2,2,1);
plot(t1, N1);
xlabel("Time [s]");
ylabel("Number of Nuclei");
title('Dt=1');

subplot(2,2,2);
plot(t5, N5);
xlabel("Time [s]");
ylabel("Number of Nuclei");
title('Dt=5');

subplot(2,2,3);
plot(t10, N10);
xlabel("Time [s]");
ylabel("Number of Nuclei");
title('Dt=10');

subplot(2,2,4);
plot(t20, N20);
xlabel("Time [s]");
ylabel("Number of Nuclei");
title('Dt=20');


% define the function
function [N, time] = euler_decay_ode(gamma, ttot, Dt, Nini)
	arguments
		gamma (1,1) double {mustBeReal, mustBePositive}
		ttot (1,1) double {mustBeReal, mustBePositive}
		Dt (1,1) double {mustBeReal, mustBePositive} = 0.01
		Nini (1,1) double {mustBeReal, mustBePositive} = 1
	end
	% solve a basic radioactive decay ode with Euler's method.

	% calculate parameters
	steps = ttot/Dt;

	% define decay ode
	decayode = @(N) -gamma.*N;

	% for loop over decay ode
	N = zeros(steps, 1);
	time = 0:Dt:ttot;
	N(1) = Nini;
	for i = 1:steps
		N(i+1) = N(i) + decayode(N(i))*Dt;
	end
end