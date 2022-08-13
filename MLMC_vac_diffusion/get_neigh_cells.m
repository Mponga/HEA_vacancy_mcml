function [neigh_cells] = get_neigh_cells(p0, xdim, ydim, zdim, lx, ly, lz)

    %% Return neighbor cells index number
    
    neigh_cells = zeros(27,3);
    
    [xpos, ypos, zpos ] = return_cell(p0, xdim, ydim, zdim);

    m = 1;
    for i=-1:1
        for j=-1:1
            for k=-1:1
                a = xpos+i;
                if a==0 
                    a = lx;
                elseif a>lx
                    a = 1;
                end
                b = ypos+j;
                if b==0 
                    b = ly;
                elseif b>ly
                    b = 1;
                end
                c = zpos+k;
                if c==0 
                    c = lz;
                elseif c>lz
                    c = 1;
                end
                neigh_cells(m,:) = [a b c];
                m=m+1;
            end
        end
    end
    
end

