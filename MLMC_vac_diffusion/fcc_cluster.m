function [Xcoord,Ycoord,Zcoord,Na] = fcc_cluster(Pars)

Lx=Pars.Nx;
Ly=Pars.Ny;
Lz=Pars.Nz;
a=Pars.a0;

% Lx, Ly, Lz are the dimension of the box 
% a - lattice parameter
% output parameters: Xcoord, Ycoord, Zcoord - coordinates of atoms

c = a; %for tetragonally distorted lattices c<>a
cnt = 0;                 %for counting
xyz = ones((4*Lx)*(4*Ly)*(4*Lz),3);   %preallocating of memory for xyz array in order to increase speed of calculation

for i=0:Lx
    for j=0:Ly
        for k=0:Lz
            p0 = [a*i, a*j, c*k];
            if p0(1) < Lx && p0(2) < Ly && p0(3) < Lz
                cnt = cnt + 1;
                xyz(cnt,:) = p0;
            end
            %% Atoms in the face
            p0 = [a*(i+0.5), a*(j+0.5), c*k];
            if p0(1) < Lx && p0(2) < Ly && p0(3) < Lz
                cnt = cnt + 1;
                xyz(cnt,:) = p0;
            end
            p0 = [a*(i+0.5), a*j, c*(k+0.5)];
            if p0(1) < Lx && p0(2) < Ly && p0(3) < Lz
                cnt = cnt + 1;
                xyz(cnt,:) = p0;
            end
            p0 = [a*i, a*(j+0.5), c*(k+0.5)];
            if p0(1) < Lx && p0(2) < Ly && p0(3) < Lz
                cnt = cnt + 1;
                xyz(cnt,:) = p0;
            end
        end
    end
end
Xcoord = xyz(1:cnt,1);
Ycoord = xyz(1:cnt,2);
Zcoord = xyz(1:cnt,3);
Na = length(Xcoord);

fprintf('Atoms generation successful : Na %d \n', Na)
