function [generations] = make_random_generations(Pars)

Natoms = round(Pars.Na/Pars.Ne); %Number of atoms to generate for each element

generations = zeros(Pars.Ngen,Pars.Na);

j=1;                    %init to loop over j
sites = 1:Pars.Na;      %Make a list with all the available sites.
R = zeros(1,Pars.Na);   %Store the random indexes from 1 to Nelements

MySites = zeros(Pars.Ne,Pars.Na);

for j = 1:Pars.Ngen
    MySites(j,:) = sites;
end

%%

for k = 1:Pars.Ngen

    remaining_sites = sites;

    for j = 1:Pars.Ne
        all_sites = sites;
        r1 = randsample(remaining_sites,Natoms); % Idx for element #j
        if j==1
            atoms_2_remove = r1;
        else
            atoms_2_remove = [atoms_2_remove r1];
        end
        all_sites(atoms_2_remove) = [];
        remaining_sites = all_sites; % Remove indexes in r1 to make selection from remaining sites
        R(r1)=j;
    end
    generations(k,1:Pars.Na)= R;
end

fprintf('Random generation successful : %d generations \n', Pars.Ngen)

end

