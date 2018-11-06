function C= Convert_Axis_2_AttitudeMatrix(RA,DEC,angle)
% Connvert RA,DEC into direction vector:
axis =[cosd(DEC)* cosd(RA) cosd(DEC)* sind(RA)  sind(DEC)];

% Calculate rotation matrix

% The z axis of camera
c3= axis;

% The x axis of camera
c1x=  c3(2)/sqrt((c3(1)^2)+(c3(2)^2));
c1y= -c3(1)/sqrt((c3(1)^2)+(c3(2)^2));
c1z=0;
c1= [c1x c1y c1z];

% The y axis of camera
c2= cross(c3,c1);

%Conversion matrix
C=[c1;c2;c3];
