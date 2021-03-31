
function wp = waypoints
%{
ME528 Group Assignment Waypoint Generator

Function that generates 5 waypoints. The waypoints are randomly generated
within certain limits each time the function is called. 
The UAV will start at the first waypoint (given by wp(1,1:3)) and must pass
each subsequent waypoint in order from row 1 to row 5. The UAV must pass 
within a distance of 1 m of each waypoint. The first waypoint will always
be (0,0,100) m. 

INPUTS: none

OUTPUTS:
wp = matrix of waypoints, where wp(i,1:3) are the x,y,z coordinates of the ith
waypoint in the inertial reference frame. The z-coordinate is height,
relative to the ground. All coordinates are given in metres. 
%}

rng('shuffle');

% initialise matrix to starting point of the UAV
wp = meshgrid([0 0 100], 1:5);

% Second waypoint
wp(2,1) = 500 + rand(1)*100;
wp(2,2) = wp(2,1)*tan(pi/6 + rand(1)*0.2);

% Third waypoint
dy3 = 800 + rand(1)*200;
an3 = deg2rad(70) +  rand(1)*deg2rad(15);
dx3 =1000*cos(an3);
wp(3,1:2) = [wp(2,1)-dx3, wp(2,2)+dy3];

% Fourth waypoint
dx4 = 1800 + rand(1)*500;
an4 = deg2rad(10) +  rand(1)*deg2rad(15);
dy4 = 1000*sin(an4);
wp(4,1:2) = [wp(3,1)+dx4, wp(3,2)+dy4];

% Fifth and last waypoint
dx5 = 100 - rand(1)*20;
an5 = deg2rad(15) - rand(1)*deg2rad(30);
dy5 = (1000 + rand(1)*500)*cos(an5);
wp(5,1:2) = [wp(4,1)+dx5, wp(4,2)+dy5];

return

