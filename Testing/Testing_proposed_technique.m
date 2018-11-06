% This script is to test the proposed star recognition algorithm

clear all;

clc;

close all;

clear classes;


% Initialize camera parameter

% SI STS parameters

FOV = 15;

img_height = 1024;

img_width = 1024;

pixel_size = 13.3/1024;

f = (img_height)*pixel_size /2/ tand(FOV/2);

% g=4;
%
% pr=ceil(384/g);

img_dimension=[img_width img_height];

%Load SPD

% pattern_database_file='pattern_catalogue_proposed_v4_v384.txt';

%Load SC

catalog_path='..\simulate\SKY2000_Magnitude6_doublestars_0.12.txt';

% Initialize camera attitude

int=2; %2.5: 10368, 3.75: 4608, 5:2592, 7.5: 1152

RAini=1;%from [0, 360];

RAend=359;

RA=[RAini:int:RAend]'; %row-wise

DECini=-89; %from(-90,90);

DECend=89;

DEC=[DECini:int:DECend]; %Column-wise

angle=45;

var_test = 0;

%feature_vector tolerance

% tol=1;

% fn_threshold=0.7;


%Simulator parameters

cent_variance=0;

no_ran_star=0; % Uncomment lines 268 - 283 if you add this

SNR=10; % SNR = 10 : ideal.

background_noise=0.0; % background_noise = 0.0 : ideal

PSF_set=0;          % Change this setting to (1,2 or 3) if you want to calculate the centroid. PSF_set = 0 means exact co-ordinate of the star given with no magnitude consideration.

no_miss = 0;

%centroider parameters

thresh = 0.3;

deviate_limit =0.0;

% Read star coordinates in Earth reference frame from star catalog

[SKYMAP_No,star_RA,star_DEC,star_MAG]= textread(catalog_path,'%d %f %f %f');


% Si is coordiante of star in Earth reference frame,

% the 3 column are X, Y,and Z

Si = [cosd(star_DEC).*cosd(star_RA) cosd(star_DEC).*sind(star_RA)  sind(star_DEC)];

% for i=1: length(star_RA)

%     % Convert RA, and DEC into ECI unit vector

%     ECI_vector =[cosd(star_DEC(i))* cosd(star_RA(i)) cosd(star_DEC(i))* sind(star_RA(i))  sind(star_DEC(i))];

%     Si(i,:)= ECI_vector;

% end

catalog=struct('SKYMAP_No',SKYMAP_No,'star_RA',star_RA,'star_DEC',star_DEC,'star_MAG',star_MAG,'Si',Si);

% Testing   variables
acc = 0;
match_acc=0;
count_acc = 0;
time_elapsed_algo = [];
time_elapsed_SPD = [];
Star_ref_empty_cases = 0;
Star_ref_empty_cases_list = [];
not_successful_RA_empty = [];
not_successful_DEC_empty = [];
multiple_star_id = 0;
failed_images_number_of_stars = [];
false_match_RA = [];
false_match_DEC = [];
avg_number_star_return = [];
missing_stars_list = [];
x = 1;
bin_size = 512;
sp_pho_top_50_2_count = [];
sp_pho_top_50_max = [];
sp_pho_id = [];
match_check = [];

% starID=zeros(size(RA,1),size(DEC,2));

% star_q_fov=zeros(size(RA,1),size(DEC,2));

% star_q_ad=zeros(size(RA,1),size(DEC,2));

% star_num_returned_match=zeros(size(RA,1),size(DEC,2));

% starnum_ad_frm_sra=zeros(size(RA,1),size(DEC,2));

% starnum_fov_frm_ad=zeros(size(RA,1),size(DEC,2));

%% THE LOOP FOR GENERATING IMAGES


for ii = 1:size(RA,1)
    
    for jj = 1:size(DEC,2)
        
        id=[];
        
        fn=[];
        
        temp=[];
        
%                 RA = 3.0;
%                 DEC = -53.0;
        % Generate sky image at predetermined attitude
        
        Reci2body= Convert_Axis_2_AttitudeMatrix(RA(ii),DEC(jj),angle);
        
        [star_matrix, I]= Plot_sky_images_mine(Reci2body, FOV, img_height, img_width, pixel_size,cent_variance, no_ran_star, SNR, background_noise, PSF_set,catalog);
        
        [row_star_matrix,col_star_matrix] = size(star_matrix);
        
        for x = 1:row_star_matrix
            id = [id star_matrix(x,1)];
        end
        
        id_centroid  =[];
        
        for x = 1 : size(star_matrix,1)
            id_centroid = [id_centroid ; star_matrix(x,1) star_matrix(x,11) star_matrix(x,12)];
        end
        
        len_id = length(id);
        
        %         if (len_id > 3)
        
        id = [];
        
        %             imshow(I);
        
        % --- This is for exact co-ordinates of the star with PSF_setting = 0 ----
        
        mag_list_sorted = [];
        mag_list = [];
        new_mag_list = [];
        centroid = [];
        centroid_old = [];
        Sc = [];
        match = 0;
        for i = 1 : len_id
            mag_list = [mag_list star_matrix(i,4)];
        end
        
        mag_list_sorted = sort(mag_list);
        
        % Missing stars due to magnitude uncertainty.
        for i = 1 : (len_id - no_miss)
            new_mag_list = [new_mag_list mag_list_sorted(i)];
        end
        
        A = (1 : len_id)' ;
        for i = 1 : len_id
            if(ismember(star_matrix(i,4),new_mag_list) == 1)
                centroid_old(i,1) = star_matrix(i,11);
                centroid_old(i,2) = star_matrix(i,12);
                Sc(i,1) = star_matrix(i,8);
                Sc(i,2) = star_matrix(i,9);
                Sc(i,3) = star_matrix(i,10);
            end
        end
        
        count_cent = 1;
        for i = 1 : length(centroid_old)
            if(centroid_old(i,1) ~= 0 && centroid_old(i,2) ~= 0)
                centroid(count_cent,1) =  centroid_old(i,1);
                centroid(count_cent,2)  = centroid_old(i,2);
                count_cent = count_cent + 1;
            end
        end
        
        %     if(length(centroid) <= 2)
        %         continue;
        %     end
        init_dist = []; % Dont mind this. Only for random numbers generation
        d=zeros(size(centroid,1),1);
        for i=1:length(centroid)
            temp= centroid(i,:)-[img_height/2 img_width/2];
            d(i)= sum(temp.*temp);
        end
        k_id= find(d==min(d));
        k_id=k_id(1);
        for i = 1 : size(centroid,1)
            dim_x = abs(centroid(i,1) - (img_height/2));
            dim_y = abs(centroid(i,2) - (img_width/2));
            
            dist = sqrt((dim_x).^2 + (dim_y).^2);
            init_dist = [init_dist ; dist];
        end
        
        % Construct the new id
        
        for x = 1 : size(centroid,1)
            for y = 1 : size(id_centroid)
                if(id_centroid(y,2) == centroid(x,1) && id_centroid(y,3) == centroid(x,2))
                    id = [id id_centroid(y,1)];
                end
            end
        end
%         id
        
        distance = [];
        
        % Adding false stars to the image - Just append the random centroid
%                 count_rand = 0; % Counting the number of random coords
%         
%                for j = 1 : 100
%                     x_coord_rand = (img_width - 0).*rand(1,1) + 0;
%                     y_coord_rand = (img_height - 0).*rand(1,1) + 0;
%                     dim_x = abs(x_coord_rand - (img_width/2));
%                     dim_y = abs(y_coord_rand - (img_width/2));
%                     d_rand = sqrt((dim_x).^2 + (dim_y).^2);
%                     if(d_rand > min(init_dist))
%                         count_rand = count_rand + 1;
%                         centroid(size(centroid,1) + 1,:) = [x_coord_rand y_coord_rand];
%                     end
%                     if(count_rand == no_ran_star)
%                         break;
%                     end
%                 end
                
        % Applying positional deviation.
        
        for i = 1 : size(centroid,1)
            distance(i) = sqrt((centroid(i,1) - centroid(k_id,1)) ^ 2 + (centroid(i,2) - centroid(k_id,2)) ^ 2);
        end
        
        ind_dist_max = find(distance == max(distance));
        
        for i = 1 : size(centroid,1)
            if(i ~= k_id)
                centroid(i,1) = ((-1).^ i) * ((distance(i) / max(distance)) * (deviate_limit)) +  (centroid(i,1));
                centroid(i,2) = ((-1).^ i) *((distance(i) / max(distance)) * (deviate_limit)) + (centroid(i,2));
            end
            if(i == k_id)
                centroid(i,1) = centroid(i,1);
                centroid(i,2) = centroid(i,2);
            end
        end
        
        tic
        
        % SPD Read:
        
        SPD_top_4_dist = dlmread('..\SPD\SPD_top_4_dist_list_Mv_6.txt');
        SPD_vect_patt_dist_1 = dlmread('..\SPD\SPD_vect_patt_Mv_6_dist_1.txt');
        SPD_vect_patt_dist_2 = dlmread('..\SPD\SPD_vect_patt_Mv_6_dist_2.txt');
        SPD_vect_patt_dist_3 = dlmread('..\SPD\SPD_vect_patt_Mv_6_dist_3.txt');
        SPD_vect_patt_dist_4 = dlmread('..\SPD\SPD_vect_patt_Mv_6_dist_4.txt');
        
        toc
        
        time_elapsed_SPD = [time_elapsed_SPD ; toc];
        
        tic
        % Co-ordinates of the star nearest to the center of the image.
        Star_ref = []; % Reference star found in the image
        dist_list = []; % Neighboring stars distance list
        sorted_dist_list = []; % Sorted distance list from nearest to the farthest neighboring star from the reference star
        top_4_dist_list = []; % Nearest 4 distances to the reference star in the image
        dist_cord_list = []; % List containing the distances and their corresponding co-ordiantes
        sorted_cord_list = []; % Sorted co-ordinates list
        ind_vector_list = []; % Individual vectors
        ind_vector_norm_list = []; % Normalized individual vectors
        bound_vector_list = []; % Boundary vectors
        bound_vector_norm_list = []; % Normalized boundary vectors
        bound_vector_new_list = []; % New frame boundary vectors list
        bound_vector_new_norm_list = []; % Normalized new boundary vectors list
        vect_pattern = []; % Stores the vector pattern extracted from the image
        inter_mat = []; % Stores the intermediate matching data > Format for the store is : (id) (ind_dist_image) (ind_dist_SPD)
        dist_tol = 4; % Distance matching tolernace (Keep around 8 for the real images)
        vect_tol = 8; % Vector pattern match tolerance (Keep around 6 for the real images)
        voting_list = []; % Exclusive voting list containing all the votes of the vector pattern matches
        
        for i = 1 : size(centroid,1)
            dim_x = abs(centroid(i,1) - (img_height/2));
            dim_y = abs(centroid(i,2) - (img_width/2));
            
            dist = sqrt((dim_x).^2 + (dim_y).^2);
            dist_list = [dist_list ; dist];
        end
        
        ind = find(dist_list == min(dist_list));
        
        center_x = centroid(ind,1);
        center_y = centroid(ind,2);
        
        dist_list = [];
        
        for i = 1 : size(centroid,1)
            if(i ~= ind)
                dim_x = abs(centroid(i,1) - center_x);
                dim_y = abs(centroid(i,2) - center_y);
                dist = sqrt((dim_x).^2 + (dim_y).^2);
                if(dist <= 510)
                    dist_cord_list = [dist_cord_list; dist centroid(i,1) centroid(i,2)];
                    dist_list = [dist_list ; dist];
                end
            end
        end
        
        sorted_dist_list = sort(dist_list);
        
        % Put the condition that the number of stars in the image must be greater than or equal to (Put this w.r.t to the dist_list)
        
        if(length(sorted_dist_list) >= 5)
            
            % Getting the nearest 4 distances to the reference star.
            
            for j = 1 : 4
                top_4_dist_list = [top_4_dist_list ; sorted_dist_list(j)];
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
                vect = [(sorted_cord_list(j,1) - center_x) (sorted_cord_list(j,2) - center_y)];
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
            
            % Converting the vectors into the local frame.
            % 4 pattern vectors have to be formed starting the traveling from different points
            
            
            for ind_dist = 1 : 4
                if(ind_dist ~= 1)
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
                vect_pattern(ind_dist,1) = bound_vector_new_list(1,1); % First element remains the same
                vect_pattern(ind_dist,2) = bound_vector_new_list(1,2); % First element remains the same                
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
                    
                    vect_pattern(ind_dist,vect_patt_ind) = vect_pattern_x;
                    vect_pattern(ind_dist,vect_patt_ind + 1) = vect_pattern_y;
                    vect_patt_ind = vect_patt_ind + 2;
                end
                bound_vector_new_list = [];
                bound_vector_new_norm_list = [];
            end
            
            
            % Matching begins from here.
            
            % Firstly, checking which distances match
            % Only check the first four distances --- Put all the matches into the inter matrix
            
            for ind_SPD = 1 : size(SPD_top_4_dist,1)
                for j = 1 : size(top_4_dist_list,1)
                    for k = 1 : 4 % Fixed size of the SPD --- Change this if you change the number of nearest taken into consideration
                        if((top_4_dist_list(j) - dist_tol) <= SPD_top_4_dist(ind_SPD,k) && SPD_top_4_dist(ind_SPD,k) <= (top_4_dist_list(j) + dist_tol))
                            inter_mat = [inter_mat ; ind_SPD j k];
                        end
                    end
                end
            end
            
            % Checking the vector pattern match and voting:
            for j = 1 : size(inter_mat,1)
                vote = 0;
                % Image vector pattern
                img_vect_patt = vect_pattern(inter_mat(j,2),:);
                img_vect_patt(img_vect_patt == 0) = [];
                %SPD vector pattern
                if(inter_mat(j,3) == 1)
                    spd_vect_patt = SPD_vect_patt_dist_1(inter_mat(j,1),:);
                end
                if(inter_mat(j,3) == 2)
                    spd_vect_patt = SPD_vect_patt_dist_2(inter_mat(j,1),:);
                end
                if(inter_mat(j,3) == 3)
                    spd_vect_patt = SPD_vect_patt_dist_3(inter_mat(j,1),:);
                end
                if(inter_mat(j,3) == 4)
                    spd_vect_patt = SPD_vect_patt_dist_4(inter_mat(j,1),:);
                end
                
                spd_vect_patt(spd_vect_patt == 0) = [];
                
                % Counting the number of matched vectors
                ind_SPD = 1;
                for k = 1 : 2 : length(img_vect_patt) - 1
                    for l = ind_SPD : 2 : length(spd_vect_patt) - 1
                        if((img_vect_patt(k) - vect_tol) <= spd_vect_patt(l) && spd_vect_patt(l) <= (img_vect_patt(k) + vect_tol) && (img_vect_patt(k+1) - vect_tol) <= spd_vect_patt(l+1)&& spd_vect_patt(l+1) <= (img_vect_patt(k+1) + vect_tol))
                            vote = vote + 1;
                            %                     ind_SPD = ind_SPD + 1;
                            break;
                        end
                    end
%                     if(vote == 0)
%                         break;
%                     end
                end
                
                % Put the vote in the corresponding inter_mat
                voting_list = [voting_list ; vote];
                inter_mat(j,4) = vote;
            end
            
            ind_max_vote = find(voting_list == max(voting_list));
            
            Star_ref = inter_mat(ind_max_vote,1);
            
            toc
            
            time_elapsed_algo = [time_elapsed_algo ; toc];
            
            %             if(ismember(Star_ref,id) == 1)
            if(Star_ref == id(ind))
                match_acc = match_acc + 1;
            else
                false_match_RA = [false_match_RA ; RA(ii)];
                false_match_DEC = [false_match_DEC ; DEC(jj)];
            end
            
            count_acc = count_acc + 1;
        end
        if(count_acc == 1000)
            return;
        end
        pause(0.2);
        
    end
end
