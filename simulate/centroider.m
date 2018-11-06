%% Centroider
% This function purposed is to determine centroid of stars in images taken
% by CMOS sensor.
% Input : grayscale image I(8 bit resolution) 
% Output: Sc star coordinate (x,y) in image plane
%         magnitude: brightness of stars
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Sc centroid magnitude] = centroider(I, FOV, img_height, img_width, pixel_size)

sizeI= size(I);

% Define gaussian filter;        
%H= fspecial('gaussian');
% Apply gaussian filter into image
%G=imfilter(I,H);

%Apply median filtering
M = medfilt2(I);

% Apply level thresholding grayscale image to BW images, level= 0.85;
level = 0.3; % Use optimum threshold
T= im2bw(M, level);

% Create structuring element:
se= strel('disk',3);
morpho_close= imclose(T,se);
%imview(morpho_close);

% Clear unused images
clear H;
clear G;
clear T;


% Label connected component , 
% This function is to label every star with an index integer number
% num : number of stars in images
% L : stars images 
[L,num]=bwlabel(morpho_close,8);

%Calculate centroid of stars                   -- Without sub-pixel
                                                                    %accuracy
% area =[];
% centroid = [];
% for i=1:num
%     [r c]= find(L==i);
%     % Calculate star area, which represent star magnitude
%     area = [area ;length(r)];
%     % Calculate star centroids
%     centroid = [centroid; (sum(r)/length(r)) (sum(c)/length(r))];
% end
%     magnitude= area;
    
    
    
 %Calculate centroid of stars       --- With sub-pixel accuracy. Apply
                                                                %weighted average method
                                                                
area = [];
centroid = [];
% matrix=[];
for i=1:num
    [yy, xx] = find(L==i);
    ind=sub2ind(size(I),yy,xx);
    % Calculate star area, which represent star magnitude
    area = [area ; length(xx)];
    maxint=max(I(ind));
    meanint= mean(I(ind));
    % Calculate star centroids
    centroid = [centroid ; sum(xx.*I(ind))/sum(I(ind)) sum(yy.*I(ind))/sum(I(ind))];
%     distance=sqrt(sum(centroid.^2,2));
%     matrix=[matrix; centroid, distance, area, maxint,meanint]; 
end

magnitude = area;    
    
%   test= zeros(1024,1024);
%   for i=1: num
%       test(round(centroid(i,1))-1:round(centroid(i,1))+1, round(centroid(i,2))-1:round(centroid(i,2))+1 )=1;
%   end
%   imview(test)
%   imview(morpho_close)

%% Transform stars coordinate into star vectors Sc in Camera frame

f= (img_height)*pixel_size /2 / tand(FOV/2);
Sc      =[];
vector  =[];
for i=1:num
    x=(centroid(i,1)-img_height/2)*pixel_size;
    y=(centroid(i,2)-img_width/2) *pixel_size;
    unit= sqrt(x^2 + y^2 + f^2);
    vector = [x/unit y/unit f/unit];
    Sc= [Sc; vector ];
end

% Sc = [Sc centroid magnitude];

