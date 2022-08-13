function [xpos, ypos, zpos] = return_cell(p0,xdim,ydim,zdim);

    xpos = -1;
    ypos = -1;
    zpos = -1;
    
    xlast = length(xdim);
    ylast = length(ydim);
    zlast = length(zdim);
    
    if p0(1) == xdim(xlast)
        xpos = xlast-1;
    end
    
    if p0(2) == ydim(ylast)
        ypos = ylast-1;
    end
    
    if p0(3) == zdim(zlast)
        zpos = zlast-1;
    end
    
    for i=1:length(xdim)-1
        if p0(1) >= xdim(i) && p0(1) < xdim(i+1)
            xpos = i; 
            i=length(xdim);
        end
    end
    
    for i=1:length(ydim)-1
        if p0(2) >= ydim(i) && p0(2) < ydim(i+1)
            ypos = i; 
            i=length(ydim)-1;
        end
    end
    
    for i=1:length(zdim)-1
        if p0(3) >= zdim(i) && p0(3) < zdim(i+1)
            zpos = i; 
            i=length(zdim)-1;
        end
    end

    if xpos == -1 || ypos == -1 || zpos == -1 
    fprintf('ERROR: Could not find a match for cell position %f %f %f\n', p0(1), p0(2), p0(3));
    end
end

