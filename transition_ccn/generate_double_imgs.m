function [imgs] = generate_double_imgs(hea_vac, hea_vac_2, Pars)

newimgs = zeros(Pars.pixel_Nx,Pars.pixel_Nx,2,length(hea_vac(:,1)));

% first_image = generate_imgs(hea_vac,Pars);
% 
% second_image = generate_imgs(hea_vac_2,Pars);

first_image = generate_large_imgs(hea_vac,Pars);

second_image = generate_large_imgs(hea_vac_2,Pars);

newimgs(:,:,1,:) = first_image;
newimgs(:,:,2,:) = second_image;

imgs = newimgs;

end

