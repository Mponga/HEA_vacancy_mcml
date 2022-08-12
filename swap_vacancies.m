function [flipval, swap, kt] = swap_vacancies(Pars,vacnet1,trans_net1,mygen,neighbours,Dvac,min_vac,vfe_val)

% Need to get the neighbor list and flip the idx there. Then generate imgs
% Get the nearest neighbors around the vacancies
nu0 = 1; % Hz: Debye frequency
NN1 = neighbours(Pars.vac_id,1:12);
fvalues = zeros(Pars.vac,12);  %Store the formation values of the possible sites
ef = zeros(Pars.vac,12);  %Store the forward  values of the possible sites
eb = zeros(Pars.vac,12);  %Store the backward values of the possible sites

cur_vac = [mygen(neighbours(Pars.vac_id,:)) Pars.vac_type];

% Generate the current images for the transition network
current_img = generate_large_imgs(cur_vac,Pars); 

transition_imgs(:,:,1,:) = current_img;

for j=1:12
    % This is the j-th site I am going to exchange with the vacancy
    jsite = NN1(:,j);
    
    % Get the nearest neighbors of the j-sites
    NN1j = neighbours(jsite,:);
    
    idx_2_flip = zeros(Pars.vac,1);
    
    % % Get the idx of the vacancy atom in the j-th neighbor
    for i=1:Pars.vac
        idx_2_flip(i,1) = find(NN1j(i,:) == Pars.vac_id(i));
        NN1j(i,idx_2_flip(i)) = jsite(i);
    end
    
    % Make the array, and generate the new images
    new_vac = [mygen(neighbours(NN1j)) Pars.vac_type];
    
    new_img = generate_large_imgs(new_vac,Pars); %generate_imgs(new_vac,Pars);  %Generate nvac new images only.
    
    fvalues(:,j) = (Dvac)*predict(vacnet1,new_img)+min_vac; %rescale values.

    transition_imgs(:,:,2,:) = new_img;

    tvalues = predict(trans_net1,transition_imgs);

    ef(:,j) = tvalues(:,1);
    eb(:,j) = tvalues(:,2);

end

(ef*Pars.Det)+Pars.mineb;
(eb*Pars.Det)+Pars.mineb;

%% Check change in the energy w.r.t to previous vacancy

kij  = zeros(Pars.vac,12);
prob = zeros(Pars.vac,12);
lower_energy_flag = zeros(Pars.vac,1);

for i=1:Pars.vac
    vfi = vfe_val(i);
    for j = 1:12
        vfj = fvalues(i,j);
        
        kij(i,j) = nu0*exp(-(ef(i,j)+(vfj))/Pars.kbT);

        if vfi >= vfj
            prob(i,j) = 1;
            lower_energy_flag(i,1) = 1; %if this happens, then we can proceed to move the vacancies normally without worrying too much
        else
            prob(i,j) = exp( -(vfj-vfi)/Pars.kbT );
            rho1 = rand(1);
            if prob(i,j) > rho1
                prob(i,j) = rho1;
            end
        end
        
    end
end

% According to the papers but seems to overestimate the diffusivity 
kt = sum(kij,2);

%%
flipval = zeros(Pars.vac,1);
flippos = zeros(Pars.vac,1);

swap    = zeros(Pars.vac,2);

for j=1:Pars.vac
       
    if lower_energy_flag(j) == 1
        
        [myval, mypos]=min(fvalues(j,:));
        rho1 = rand(1);
        if rho1 < 0.05 %0.05 %Add some randomness here to allow for some randon jump even when the system is supposed to minimize energy. 
            mypos = randi([1,12],1);
        end
        flipval(j,1) = myval;
        flippos(j,1) = mypos;
        swap(j,1) = NN1j(j,flippos(j,1)); %-> Vacancy to migrate here. Need to update Pars.vac_id and in mygen
        swap(j,2) = Pars.vac_id(j); %-> Vacancy to migrate here. Need to update Pars.vac_id and in mygen
    else
        
        [myval, mypos]=max(prob(j,:));
        flipval(j,1) = fvalues(j,mypos);
        flippos(j,1) = mypos;
        swap(j,1) = NN1j(j,flippos(j,1)); %-> Vacancy to migrate here. Need to update Pars.vac_id and in mygen
        swap(j,2) = Pars.vac_id(j); %-> Vacancy to migrate here. Need to update Pars.vac_id and in mygen
    end
    
end 


end

