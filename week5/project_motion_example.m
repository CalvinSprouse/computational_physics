%Projectile Motion Example

%{ 
This program calculate the time of flight of an object undergoing
projectile motion. It also plots the trajectory.
%}

%Author: Darci Snowden
%Date: January 26th, 2024

%Version 1: Example version for students to edit

%Define the acceleration of gravity as a global variable
global g;
g=1.6;  %acceleration of gravity on the Moon, m/s^2

%Ask the user for their input
angle_launch=input('What is the launch angle in degrees? ');
initial_v=input('What is the initial velocity in m/s? ' );
initial_h=input('What is the initial height in m? ' );

%Determine the impact 
tfin = findImpactTime(initial_h,initial_v,angle_launch);

%Make a plot of the trajectory
drawProjMotPlot(initial_h,initial_v,angle_launch,tfin);

% print the range of the projectile
proj_range = findRange(initial_v, angle_launch, tfin);
fprintf("The projectile travelled %.2f meters.\n", proj_range);


function drawProjMotPlot(y0,v0,launchAng,tfin)
    global g;
    
    %Define functions for graph
    xpos=@(t) v0*cosd(launchAng)*t;
    ypos=@(t) y0+v0*sind(launchAng)*t-1/2*g*t^2;

    %Plot the function defining the x position vs. the y position
    fplot(xpos,ypos,[0 tfin],'r','LineWidth',2.);
    xlabel('X Position (m)');
    ylabel('Y Position (m)');
    title('Projectile Motion Example');
    set(gca,'FontSize',14);
    
end

function t_impact = findImpactTime(y0,v0,launchAng)

    global g;

    %Find the initial velocity
    vy0=v0*sind(launchAng);
 
    %Use the quadratic formula to find the impact time
    t_impact=(-vy0-sqrt(vy0^2+2*g*y0))/(-g);
    
end

function range = findRange(v0, launchAng, tfin)
	% return the range of a projectile
	% v0 : the initial speed [m/s]
	% launchAng : the initial angle [degrees]
	% tfin : the airtime of the projectile [s]

	% calculate the x velocity
	xvel = v0 * cosd(launchAng);

	% calculate the distance traveled
	range = xvel * tfin;
end