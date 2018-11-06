% This script is to generate the catalogue for the proposed idea to prevent


clear all;
clc;
%% Initialize camera parameter
img_height = 1024;
img_width = 1024;
pixelsize= 13.3/1024;
FOV =15;
angle = 45;
test = 0;
% Max_N = 15;

file_path='..\simulate\SKY2000_Magnitude6_doublestars_0.12.txt';
[SKYMAP_No,star_RA,star_DEC,star_MAG]= textread(file_path,'%d %f %f %f');

% Si is coordiante of star in Earth reference frame,
% the 3 column are X, Y,and Z
no_stars=length(star_RA);
Si = zeros(no_stars, 3);
catalog = zeros(no_stars, 4);

% Database prepared
vect_pattern = []; % Final added vector pattern
top_4_dist = []; % Nearest 4 distances to the star id

for i=1: no_stars
    % Convert RA, and DEC into ECI unit vector
    %     i = 31;
    ECI_vector =[cosd(star_DEC(i))* cosd(star_RA(i)) cosd(star_DEC(i))* sind(star_RA(i))  sind(star_DEC(i))];
    Si(i,:)= ECI_vector;
end

for i=1: no_stars
%     i = 31; % Initial testing purposes.
    
    %Information
    RA= star_RA(i);
    DEC= star_DEC(i);
    dist_list = [];  %Only the distance information
    dist_cord_list = []; % Distance and their corresponding co-ordinates list
    sorted_dist_list = []; % Sorted distance information
    sorted_cord_list = []; % Sorted co-ordinates according to the distances from the center star.
    ind_vector_list = []; % Individual vector list
    ind_vector_norm_list = []; % Individual normalized vector list
    bound_vector_list = []; % Boundary vector list
    bound_vector_norm_list = []; % Normalized boundary vector list
    bound_vector_new_list = []; % New frame for the boundary vector list
    bound_vector_new_norm_list = []; % New frame for the normalized boundary vector list
    
    C= Convert_Axis_2_AttitudeMatrix(RA,DEC,angle);
    [R_camera_to_earth,star_matrix]= Find_neighbor_star_half_FOV(C,FOV, img_height, img_width, pixelsize);
    
    for j = 1 : size(star_matrix,1)
        
        if(star_matrix(j,1) ~= i)
            x_cord = star_matrix(j,11);
            y_cord = star_matrix(j,12);
            
            % x_cord corresponds to the row number of the star in the image.
            % y_cord corresponds to the column number of the star in the image.
            
            dim_x = abs(x_cord - (img_height/2));
            dim_y = abs(y_cord - (img_width/2));
            
            
            % Arranging the neighboring stars in the order of distances from
            % the center star
            
            dist = sqrt((dim_x).^2 + (dim_y).^2);
            dist_cord_list = [dist_cord_list; dist x_cord y_cord];
            dist_list = [dist_list ; dist];
        end
    end
    
    sorted_dist_list = sort(dist_list);
    
    % Getting the top four distances for each star entry
    
    for j = 1 : 4
        top_4_dist(i,j) = sorted_dist_list(j);
    end
    
    % Arranging the coord list according to their distance order.
    
    for j = 1 : size(sorted_dist_list,1)
        for k = 1 : size(sorted_dist_list,1)
            if(dist_cord_list(k,1) == sorted_dist_list(j))
                sorted_cord_list = [sorted_cord_list ; dist_cord_list(k,2) dist_cord_list(k,3)];
                break;
            end
        end
    end
    
    % Getting the individual vector list and the boundary vector list of the stars.
    
    for j = 1 : size(sorted_cord_list,1)
        vect = [(sorted_cord_list(j,1) - (img_width/2)) (sorted_cord_list(j,2) - (img_height/2))];
        ind_vector_list = [ind_vector_list ; vect];
        vect = vect / norm(vect);
        ind_vector_norm_list = [ind_vector_norm_list ; vect];
    end
    
    for j = 1 :  size(sorted_cord_list,1) - 1
        bound_vector = [(ind_vector_list(j+1,1) - ind_vector_list(j,1)) (ind_vector_list(j+1,2) - ind_vector_list(j,2))];
        bound_vector_list = [bound_vector_list ; bound_vector];
        bound_vector = bound_vector / norm(bound_vector);
        bound_vector_norm_list = [bound_vector_norm_list ; bound_vector];
    end
    
    % Converting the vectors into the local frame and starting from one point
    % at a time.
    
    ind_dist = 4 %% CHANGE  this for getting the SPD for the corresponding travel point start (nearest 1st star ; nearest 2nd star ; nearest 3rd star ; nearest 4th star)
    
    for j = 1 : ind_dist - 1
        ind_vector_list(1,:) = [];
        ind_vector_norm_list(1,:) = [];
        bound_vector_list(1,:) = [];
        bound_vector_norm_list(1,:) = [];
    end
    
    for j = 1 : size(bound_vector_list,1)
        bound_vector_new_x = bound_vector_list(j,:) * transpose(ind_vector_norm_list(j,:));
        bound_vector_new_y = bound_vector_list(j,:) * transpose(ind_vector_norm_list(j+1,:));
        bound_vector_new = [bound_vector_new_x  bound_vector_new_y];
        bound_vector_new_list = [bound_vector_new_list;  bound_vector_new_x bound_vector_new_y];
        bound_vector_new = bound_vector_new / norm(bound_vector_new);
        bound_vector_new_norm_list = [bound_vector_new_norm_list ; bound_vector_new];
    end
    
    % Generating the vector pattern as we travel along the path.
    vect_pattern(i,1) = bound_vector_new_list(1,1); % First element remains the same
    vect_pattern(i,2) = bound_vector_new_list(1,2); % First element remains the same
    vect_patt_ind = 3;
    
    for j = 2 : size(bound_vector_new_list,1)
        vect_pattern_x = 0;
        vect_pattern_y = 0;
        
        % Normal part
        vect_pattern_x = vect_pattern_x + bound_vector_new_list(1,1);
        vect_pattern_y = vect_pattern_y + bound_vector_new_list(j,2);
        
        %Residue part of the x_pattern
        for k = j : -1 : 2
            vect_pattern_x = vect_pattern_x + (bound_vector_list(k,:) * transpose(ind_vector_norm_list(1,:)));
        end
        
        %Residue part of the y_pattern
        for k = j-1 : -1 : 1
            vect_pattern_y = vect_pattern_y + (bound_vector_list(k,:) * transpose(ind_vector_norm_list(j+1,:)));
        end
        
        vect_pattern(i,vect_patt_ind) = vect_pattern_x;
        vect_pattern(i,vect_patt_ind + 1) = vect_pattern_y;

        vect_patt_ind = vect_patt_ind + 2;
        
    end    
end  
