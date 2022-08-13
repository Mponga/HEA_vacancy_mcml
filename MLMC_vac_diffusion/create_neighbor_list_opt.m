function [NN1,NumNeigh] = create_neighbor_list_opt(atoms, Pars, Cells)

Rcut = Pars.R0;
Nx = Pars.Nx;
Ny = Pars.Ny;
Nz = Pars.Nz;

xdim = Cells.xdim;
ydim = Cells.ydim;
zdim = Cells.zdim;
lx   = Cells.lx;
ly   = Cells.ly;
lz   = Cells.lz;
CellIdx = Cells.CellIdx;
CellArraySize = Cells.CellArraySize;

NNZ = 54;

Na = length(atoms(:,1));
NN1 = zeros(Na,NNZ);         %store 54 values for NN, increase if needed
NumNeigh = zeros(Na,1);


if Na < 1000
    [NN1,NumNeigh] = create_neighbor_list(atoms,Pars);
else
    
    parpool(4)
    parfor i=1:Na
        
        if rem(i/Na*100,10) == 0
            fprintf('Completed %d of 100 percent\n',i/Na*100);
        end
        
        v = zeros(1, NNZ);
        distance = zeros(NNZ,3);
        Rsqr = zeros(NNZ,1);
        neigh_cells = zeros(27,3);
        
        p0 = [atoms(i,1) atoms(i,2) atoms(i,3)];
        
        neigh_cells = get_neigh_cells(p0, xdim, ydim, zdim, lx, ly, lz);
        linearInd = sub2ind(CellArraySize,neigh_cells(:,1),neigh_cells(:,2),neigh_cells(:,3));
        AtomsIDToCheck = nonzeros(CellIdx(linearInd,:));
        MyAtomsToCheck = length(AtomsIDToCheck);
        
        m=1;
        for j=1:MyAtomsToCheck
            if i~=AtomsIDToCheck(j)
                p1 = [atoms(AtomsIDToCheck(j),1) atoms(AtomsIDToCheck(j),2) atoms(AtomsIDToCheck(j),3)];
                [R, dx, dy, dz] = get_distance_pbc(p1,p0,Nx,Ny,Nz);
                if R < Rcut
                    v(1,m)=AtomsIDToCheck(j);
                    
                    distance(m,:) = [dx dy dz];
                    Rsqr(m,1) = R*R;
                    
                    m=m+1;
                end
            end
        end
        
        %         [A, sortedIdx] = sortrows(distance,2);
        [A, sortedIdx] = sort(Rsqr); %want to get the NN1 first.
        NN1(i,:) = v(sortedIdx);
        
        NumNeigh(i,1) = m-1;
    end
    
    p = gcp;
    delete(p)
end

end

