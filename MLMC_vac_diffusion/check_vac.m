function [v] = check_vac(x, Pars)

out =0; %Init output to error in case we want to print a msg.

atoms = x(Pars.vac_id,:);

Rcut=Pars.R0/2;
Nx=Pars.Nx;
Ny=Pars.Ny;
Nz=Pars.Nz;

v = ones(Pars.vac,1); %Assume all vacancies are available

for i=1:Pars.vac
    
    p0 = [atoms(i,1) atoms(i,2) atoms(i,3)]; %Get one vacancy
    
    for j=1:Pars.vac
        if i ~= j && v(i) ~= 0
            p1 = [atoms(j,1) atoms(j,2) atoms(j,3)];
            R = get_distance_pbc(p1,p0,Nx,Ny,Nz);
            if R <= Rcut
                v(i)=0; %Vacancy are too close, make then unavailable.
                v(j)=0; %Vacancy are too close, make then unavailable.
            end
        end
    end
end

%%
out =1; %set to success in case we want to print a msg.

end

