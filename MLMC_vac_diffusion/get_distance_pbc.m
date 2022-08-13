function [R, dx, dy, dz] = get_distance_pbc(a,b,xsize,ysize,zsize)
% Compute the distance between two vectors a and b when pbc are applied

dx = a(1)-b(1);
dy = a(2)-b(2);
dz = a(3)-b(3);

%Check PBC
if dx >  0.5*xsize
    dx = dx - xsize;
    
elseif dx <= -0.5*xsize
    dx = dx + xsize;
    
end

if dy >  0.5*ysize
    dy = dy - ysize;
    
elseif dy <= -0.5*ysize
    dy = dy + ysize;
    
end

if dz >  0.5*zsize
    dz = dz - zsize;
    
elseif dz <= -0.5*zsize
    dz = dz + zsize;
    
end


R = sqrt(dx^2+dy^2+dz^2);
end

