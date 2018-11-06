%% This script is to plot stars images on sky
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function [ R_camera_to_earth, star_matrix]= Find_neighbor_star_FOV(C,FOV, img_height, img_width, pixel_size)

f= (img_height)*pixel_size /2 / tand(FOV/2);  % Change this too if you are changing the FOV size ahead.

%% Rotation matrix from  Earth to camera reference frame
R_camera_to_earth = C;

% i1=[1 0 0]';
% i2=[0 1 0]';
% i3=[0 0 1]';
% 
% c1= R_earth_to_camera * i1;
% c2= R_earth_to_camera * i2;
% c3= R_earth_to_camera * i3;
% 
% c=[c1 c2 c3];

%R_camera_to_earth = [c1' ;c2'; c3'];

c1= R_camera_to_earth(1,:);
c2= R_camera_to_earth(2,:);
c3= R_camera_to_earth(3,:); % Z-axis of the camera attitude matrix

%% Read star coordinates in Earth reference frame from star catalog
file_path='spr_ham_spear/simulate/SKY2000_Magnitude6_doublestars_0.12.txt';
[SKYMAP_No,star_RA,star_DEC,star_MAG]= textread(file_path,'%d %f %f %f');

% Si is coordiante of star in Earth reference frame, 
% the 3 column are X, Y,and Z
Si = zeros(length(star_RA), 3);
for i=1: length(star_RA)
    % Convert RA, and DEC into ECI unit vector
    ECI_vector =[cosd(star_DEC(i))* cosd(star_RA(i)) cosd(star_DEC(i))* sind(star_RA(i))  sind(star_DEC(i))];
    Si(i,:)= ECI_vector;
end

%% Find all the neighbor stars to reference stars within the half FOV
Si_FOV=[];
star_matrix = [];% This star matrix store all stars that is within FOV of camera[SKYMAP_No  Xc  Yc  Zc]
Star_ID = [];
S11 = [];
for i=1: length(star_RA)
    if dot(Si(i,:), c3') > cosd(FOV/2)      % Check this condition: This is for stars which lie within FOV / 2
        Si_FOV= [Si_FOV ; Si(i,:)];
       Star_ID = [Star_ID ; SKYMAP_No(i)];
       star_matrix = [star_matrix ; SKYMAP_No(i) star_RA(i) star_DEC(i)  star_MAG(i)];
    end
end


%% Project star into camera frame:
Sc_FOV =[];
for i=1: size(star_matrix,1)
    Sc_FOV= [Sc_FOV; dot(Si_FOV(i,:),c1') dot(Si_FOV(i,:), c2') dot(Si_FOV(i,:), c3')];
end

x=[];
y=[];
% Project stars into image frame
for i=1: size(star_matrix, 1)
    x= [x; f* Sc_FOV(i,1)/Sc_FOV(i,3)/pixel_size  + img_height/2];
    y= [y; f* Sc_FOV(i,2)/Sc_FOV(i,3)/pixel_size + img_width/2];
end

star_matrix= [star_matrix Si_FOV Sc_FOV x y];
