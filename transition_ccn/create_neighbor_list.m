function [NN1,NumNeigh] = create_neighbor_list(atoms,Pars)

Rcut=Pars.R0;
Nx=Pars.Nx;
Ny=Pars.Ny;
Nz=Pars.Nz; 

NNZ = 54;

Na = length(atoms(:,1));
NN1 = zeros(Na,NNZ);         %store 60 values for NN, increase if needed
NumNeigh = zeros(Na,1);

% parpool(3)
% parfor i=1:Na 
for i=1:Na
      if rem(i/Na*100,10) == 0 
          fprintf('Completed %d of 100 percent\n',i/Na*100);
      end
      v = zeros(1, NNZ);
        p0 = [atoms(i,1) atoms(i,2) atoms(i,3)];
        m=1;
        for j=1:Na
            if i~=j
                p1 = [atoms(j,1) atoms(j,2) atoms(j,3)];
                R = get_distance_pbc(p1,p0,Nx,Ny,Nz);
                if R <= Rcut
                    v(1,m)=j;
                    m=m+1;
                end
            end
        end
        NN1(i,:)= v;
        NumNeigh(i,1) = m-1;
end

% p = gcp;
% delete(p)

fprintf('Neighbor list successful \n')


end

