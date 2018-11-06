% This script is to plot stars images on sky
%--------------------------------------------------------------------------

% INPUTS:
%  C: Rotation Matrix (eci2body)
%  FOV: Camera field of view in degrees
%  img_height: num of pixels in height
%  img_width: num of pixels in width
%  pixel_size:pixel pitch in mm.

% OUTPUTS:
% R_camera_to_earth: Rotation Matrix (eci2body)
% star_matrix: [SKYMAP_No(i) star_RA(i) star_DEC(i)  star_MAG(i) Xi Yi Zi Xb Yb Zb Ximg Yimg]
% I: simulated image

%--------------------------------------------------------------------------

function [ star_matrix I]= Plot_sky_images_mine(C, FOV, img_height, img_width, pixel_size, variance, no_ran_star, SNR, background_noise, PSF_set, catalog)

SKYMAP_No=catalog.SKYMAP_No;
star_RA=catalog.star_RA;
star_DEC=catalog.star_DEC;
star_MAG=catalog.star_MAG;

f= (img_height)*pixel_size /2 / tand(FOV/2);    % Change this for which ever FOV you are chosing.

% Rotation matrix from  Earth to camera reference frame
Reci2body = C; %This input is R(eci2body)

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

c1= Reci2body(1,:);% x-body in ECI coords
c2= Reci2body(2,:);% y-body in ECI coords
c3= Reci2body(3,:);% z-body in ECI coords



% Si is coordiante of star in Earth reference frame,
% the 3 column are X, Y,and Z
Si = catalog.Si;

% Sc = zeros(length(star_RA), 3);
% %Rotate star coordinates from Earth to Camera reference frame
% %Sc is star coordinate in Camera reference
% for i=1: length(star_RA)
%     vector = R_earth_to_camera * Si(i,:)'  ;
%     Sc(i,:)= vector';
% end


%% Find stars within FOV of camera
Si_FOV=[];
star_matrix = [];% This star matrix store all stars that is within FOV of camera[SKYMAP_No  Xc  Yc  Zc]
for i=1: length(star_RA)
    if dot(Si(i,:), c3') > cosd(FOV/2)     % Change this for varying size of the FOV
        Si_FOV= [Si_FOV ; Si(i,:)]; %Store the ECI coordinates of stars within FOV
        star_matrix = [star_matrix ; SKYMAP_No(i) star_RA(i) star_DEC(i)  star_MAG(i)]; %Include all the stars within FOV in star_matrix
    end
end

%% Project star into camera frame:
Sc_FOV =[];
for i=1: size(star_matrix,1)
    Sc_FOV= [Sc_FOV; dot(Si_FOV(i,:),c1') dot(Si_FOV(i,:), c2') dot(Si_FOV(i,:), c3')]; %Convert ECI star vectors to CamBody. Equivalent to C*Si(i,:)'
end

x=[];
y=[];
% Project stars into image frame
for i=1: size(star_matrix, 1)
    x= [x; f* Sc_FOV(i,1)/Sc_FOV(i,3)/pixel_size  + img_height/2];
    y= [y; f* Sc_FOV(i,2)/Sc_FOV(i,3)/pixel_size + img_width/2];
end

star_matrix= [star_matrix Si_FOV Sc_FOV x y];

%% Create image and add stars into image
star_coordinate_X= round(x);
star_coordinate_Y= round(y);

I= zeros(img_height, img_width);

% for i=1:size(star_matrix,1);
%     if (star_coordinate_X(i)>6) && (star_coordinate_Y(i)>6)
%     I(star_coordinate_X(i)-3:star_coordinate_X(i)+3, star_coordinate_Y(i)-3:star_coordinate_Y(i)+3)=1;
%     end
% end

%imshow(I);
%% Create image and add stars into image
% Add random stars into image
% Star stationed near the center (mean is width/2)
% Actually can just randomise within a limit, no need normal distribution.

ran_star_coor_X= random('norm', img_height/2, sqrt(img_height), 1,no_ran_star);
ran_star_coor_Y= random('norm', img_height/2, sqrt(img_height), 1,no_ran_star);
star_coordinate_X=[star_coordinate_X; ran_star_coor_X'];
star_coordinate_Y=[star_coordinate_Y; ran_star_coor_Y']; %Introducing random stars

if (length(x)>0)
    %Introduce centroiding error with given variance (due to abherration/refractive effects/sensor misalignement!)
    %Mean, sqrt(variance) in # of pixels.
    Rx=random('norm', 0, sqrt(variance), 1,size(star_matrix, 1));
    Ry=random('norm', 0, sqrt(variance), 1,size(star_matrix, 1));
    star_coordinate_X= round(x+Rx');
    star_coordinate_Y= round(y+Ry');
end


for i=1:size(star_matrix,1);
    PSF_im = PSF_mine(star_matrix(i,4), PSF_set); %Introduce PSF from defocusing based on star magnitude
    centre = (size(PSF_im,1)-1)/2;
    if (star_coordinate_X(i)>centre) && (star_coordinate_Y(i)>centre)
        I( star_coordinate_Y(i)-centre:star_coordinate_Y(i)+centre,star_coordinate_X(i)-centre:star_coordinate_X(i)+centre)=PSF_im; %Input PSF-ed image into I
    end
end

%  Read: http://stackoverflow.com/questions/16008228/using-imnoise-to-add-gaussian-noise-to-an-image
% and http://www.mathworks.com/matlabcentral/answers/24282-image-processing-noise
%imshow(I);

v_noise = var(I(:)) / (10^(SNR/10));
% Create noise image N
 I= imnoise(I, 'gaussian',background_noise,v_noise); %Introducing centroiding error due to SNR(mean, variance) in the range [0,1]. The mean(sky signal) can be found by using real images!
% I= imnoise(I, 'gaussian',0,0.0021);
% I=awgn(I,9);