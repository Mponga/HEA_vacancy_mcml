function [MyCells] = input_cells(atoms, Pars)

xcoord=atoms(:,1);
ycoord=atoms(:,2); 
zcoord=atoms(:,3);

% Create a system of cell to split atoms and reduce time for neighboud search

Nx = Pars.Nx;
Ny = Pars.Ny;
Nz = Pars.Nz;
Na = Pars.Na;

cellsize = 2;

lx = round(Nx/cellsize)-1; %Voxel size in x-direction
ly = round(Ny/cellsize)-1; %Voxel size in y-direction
lz = round(Nz/cellsize)-1; %Voxel size in z-direction

NumCells = lx*ly*lz;

Cells = zeros(lx,ly,lz);
CellIdx = zeros(NumCells,round(Na/NumCells));


m=1;
for i=1:lx-1
    for vacID=1:ly-1
        for k=1:lz-1
            Cells(i,vacID,k) = m;
            m=m+1;
        end
    end
end

xdim = linspace(0,Nx,lx+1);
ydim = linspace(0,Ny,ly+1);
zdim = linspace(0,Nz,lz+1);

%%
MyPos = ones(lx*ly*lz,1);

CellArraySize=size(Cells);

for i=1:Na
    p0 = [xcoord(i) ycoord(i) zcoord(i)];
    [xpos, ypos, zpos] = return_cell(p0, xdim, ydim, zdim);
    linearInd = sub2ind(CellArraySize,xpos,ypos,zpos);
    
    CellIdx(linearInd,MyPos(linearInd,1)) = i;
    MyPos(linearInd,1) = MyPos(linearInd,1)+1;
    
end
%% For a given position, compute the neighbor cells

% neigh_cells = get_neigh_cells(p0, xdim, ydim, zdim, lx, ly, lz);

% linearInd = sub2ind(CellArraySize,neigh_cells(:,1),neigh_cells(:,2),neigh_cells(:,3));
% AtomsIDToCheck = nonzeros(CellIdx(linearInd,:));
% MyAtomsToCheck = length(AtomsIDToCheck);

% scatter3(xcoord(AtomsIDToCheck),ycoord(AtomsIDToCheck),zcoord(AtomsIDToCheck),100,'filled')

%% Now assign values to structure

MyCells.cellsize = cellsize;

MyCells.lx = lx;
MyCells.ly = ly;
MyCells.lz = lz;

MyCells.NumCells = NumCells;
MyCells.Cells = Cells;
MyCells.CellIdx = CellIdx;

MyCells.xdim = xdim;
MyCells.ydim = ydim;
MyCells.zdim = zdim;

MyCells.MyPos = MyPos;
MyCells.CellIdx = CellIdx;

MyCells.CellArraySize = CellArraySize;

end

