clear all;
close all;

%This program demonstrates a direction field for a differential equation

%Author: Darci Snowden
%Date: October 15, 2020

%Define the position arrays x and y
[t,y] = meshgrid(0:.1:2,-2:.1:2);

%Define the rate of change array
dudt = -y;

%Create a quiver plot 
figure
quiver(t,y,t,dudt,'r')

%Plot streamlines that start at different points
startt = 0.1:.1:1;
starty = ones(size(startt));
streamline(t,y,t,dudt,startt,starty);
hold on

%Start points in the other quadrant
starty = -1*ones(size(startt));
streamline(t,y,t,dudt,startt,starty);
set(gca,'FontSize',14);
ylabel('y(t)','FontSize',14);
xlabel('t','FontSize',14);
ylim([-1.,1.])

% Set the axes to go through the origin
set(gca,'XAxisLocation','origin','YAxisLocation','origin')

grid on;
