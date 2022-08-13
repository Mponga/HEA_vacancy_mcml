function [chem_env, ave_chem, sigma] = analyze_chem_env(Pars,neighbours,NumNeigh, generations)

fprintf('Computing chemical environment \n')

chem_env = zeros(4,Pars.Ne,Pars.Na,Pars.Ngen); %Na atoms x Ne elements x 4 nighbor shells

ave_chem = zeros(4,Pars.Ne,Pars.Ngen);

for j=1:Pars.Ngen
    
    if rem(j,50) == 0
        fprintf('Analyzing generation %d of %d \n',j,Pars.Ngen);
    end
    
    for i=1:Pars.Na
        
        counts = compute_chem_env(i,neighbours,NumNeigh, generations(j,:), Pars.Ne);
        chem_env(1:4,1:Pars.Ne,i,j) = counts'./[12 18 42 54]';
        ave_chem(:,:,j) = ave_chem(:,:,j) + counts'./[12 18 42 54]';
        
    end
    
    ave_chem(:,:,j) = ave_chem(:,:,j)/Pars.Na; % Average value for NN1 to NN4 shell considering all atoms.
    
    % Compute standard deviation of ci for different neighbor shells.
    
    num = zeros(4,Pars.Ne);
    for i=1:Pars.Na
        A = (chem_env(:,:,i,j) - ave_chem(:,:,j)).^2;
        num = num + A;
    end
    
    sigma(:,:,j) = sqrt(num/Pars.Na); %Standard deviation for ci depending on the different neighbor shells.
end

fprintf('Chemical environment analysis successful \n')

end

