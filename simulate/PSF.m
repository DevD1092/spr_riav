%This script is for generating the point spread function of the stars in the simulated images.


function PSF_im = PSF(magnitude, settings)

if (settings==0)
    PSF_im=ones(1,1);
    return;
end

%% Settings 1:
if (settings==1)
    A=1;
    if magnitude <= 1
        centre=3;
        sigma=sqrt(5);
    elseif magnitude<=2
        centre=3;
        sigma=sqrt(4);
    elseif magnitude<=3
        centre=3;
        sigma=sqrt(3);
    elseif magnitude<=4
        centre=3;
        sigma=sqrt(2);
    elseif magnitude<=6
        centre=3;
        sigma=sqrt(2);
    end
    
    %% Settings 2:
elseif (settings==2)
    A=1;
    if magnitude <= 1
        centre=5;
        sigma=sqrt(10);
    elseif magnitude<=2
        centre=4.5;
        sigma=sqrt(8);
    elseif magnitude<=3
        centre=4;
        sigma=sqrt(7);
    elseif magnitude<=4
        centre=3.5;
        sigma=sqrt(6);
    elseif magnitude<=6
        centre=3;
        sigma=sqrt(4);
    end
elseif (settings==3)
    centre=10;
    if magnitude <= 1
        sigma=sqrt(8);
        A=1;
    elseif magnitude<=2
        sigma=sqrt(7);
        A=0.9;
    elseif magnitude<=3
        sigma=sqrt(5);
        A=0.8;
    elseif magnitude<=4
        sigma=sqrt(4);
        A=0.7;
    elseif magnitude<=5
        sigma=sqrt(3);
        A=0.6;
    elseif magnitude<=6
        sigma=sqrt(2);
        A=0.5;
        %A should be a fraction based on star brightness and visual magnitude.
        %A=10^(Mv/-2.5) x some multiplier
    end
    
end




height= 2*(centre-1)+1;
width= 2* (centre-1)+1;
PSF_im=zeros(height);

for x=1:height
    for y=1:width
        PSF_im(x,y)=  A *exp((-1)*((x-centre)^2+(y-centre)^2)/(2*sigma^2));
    end
end
% PSF_im= uint8(PSF_im.*255);
% imview(PSF_im);
