% get the coordinates for the surface of a sphere
[X, Y, Z] = sphere;

% set the radius of each sphere
r1 = 1;
r2 = .5;

% define the direction of motion for each sphere
d1 = -1;
d2 = 0;

% avoid multi-collision
has_collided = 0;

for xs = -5 : .1 : 5
    % clear and hold the figure window
    clf;
    hold on;

	% define sphere positions
	s1x = r1*(X + d1*xs);
	s1y = r1*Y;
	s1z = r1*Z;

	s2x = r2*(X + d2*xs);
	s2y = r2*Y;
	s2z = r2*Z;
    
    % use surface plots to plot both spheres
    surf(s1x, s1y, s1z, FaceColor=[0, 0, 1], EdgeColor='none')
    surf(s2x, s2y, s2z, FaceColor=[1, 0, 1], EdgeColor='none')

	% check for sphere contact
	% first calculate distance
	dist = sqrt( xs.^2 );

	% compare distance to radius
	if dist < (r1 + r2) && (has_collided == 0)
		% change motion direction of each sphere
		d2 = -1;
		d1 = 1;

		% store collision
		has_collided = 1;
	end
    
	% disable holding
	hold off;
    
    % stuff that controls the axes and view angle of the plot
    axis([-5, 5, -5, 5, -5, 5]);
    axis square;
    view([-20  20]);
   
    % update the image
    drawnow;
end
