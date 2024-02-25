%Incomplete 2nd Order Runge-Kutta 

clear all;
close all;


%Define ode 
func=@(t,y) -1.2*y+7*exp(-0.3*t);

%Define initial values and time step
y0=3;
tmin=0;
tmax=8;
h=0.5;

%Call 2nd Order Runge Kutta function to solve ODE
[time,y] = rk2(func,tmax,tmin,h,y0)


%Define exact solution

yexact=@(time) 70/9*exp(-0.3*time)-43/9*exp(-1.2*time);

%Plot the results
hold on;
plot(time,y,'.','MarkerSize',15.);
fplot(yexact,[tmin tmax],'r');
set(gca,'FontSize',14.);
xlabel('Time (s)','FontSize',14.);
ylabel('y(t)','FontSize',14.);
legend('Numerical Solution','Exact Solution');
hold off;

function [t,y] = rk2(func,tmax,tmin,h,y0)

%calc number of times to go through loop
N=(tmax-tmin)/h;

%initialize arrays for time and y values
t(1)=tmin
y(1)=y0

%for loop to fill rest of y and t arrays
for i=1:N 
   
    
    %calc k1 and k2
    k1 = func(t(i),y(i));
    k2 = func(t(i)+0.5*h,y(i)+0.5*h*k1);
    
    %update t and y arrays
    t(i+1)=t(i)+h;;
    y(i+1)=y(i)+h*k2;
end
end

% function [t,y] = rk2(%to fill)
% 
% %calc number of times to go through loop
% 
% 
% %initialize arrays for time and y values
% 
% 
% %for loop to fill rest of y and t arrays
% for i=1:N 
%    
%     
%     %calc ks
%     k1 = 
%     k2 = 
%     k3 =
%     k4 =
%     
%     %update t and y arrays
%     t(i+1)=
%     y(i+1)=
% end
end