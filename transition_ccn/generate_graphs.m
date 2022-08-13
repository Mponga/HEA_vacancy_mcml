function [adjacency] = generate_graphs(hea_vac, Pars, vacIDs, neighbours, typeID, ratoms, Lx, Ly, Lz)

tol = 0.25; 
% Make chemical tensor
Nvac = length(hea_vac(:,1));
Nx = Pars.Nx;
Ny = Pars.Ny;
Nz = Pars.Nz;
R0 = Pars.R0;
NN4 = 54;
rcut = 4*(Pars.latconst*sqrt(2)/2)+tol;
b = Pars.latconst*sqrt(2)/2;

% NN1 = b; NN2 = sqrt(2)*b; NN3=sqrt(3)*b; NN4=2*b;
% 1: Co, 2: Cr, 3: Fe, 4: Ni
% 1-1: 1, 2-2: 2, 3-3: 3, 4-4:4
% 1-2: 5, 1-3: 6, 1-4: 7
% 2-1: 5, 2-2: 2, 2-3: 8, 2-4: 9
% 3-1: 6, 3-2: 8, 3-3: 3, 3-4: 10
% 4-1: 7, 4-2: 9, 4-3: 10, 4-4:4

%Interaction values are well defined levels to make adjacency well
%contrasted matrix 
IntValues = [1 5 6 7;
             5 2 8 9;
             6 8 3 10;
             7 9 10 4];

mu = [8.126 7.746 7.318 7.573]; % Chemical potential from JAC. 
maxmu = max(mu);
minmu = min(mu);
deltamu = maxmu-minmu;
barmu = sum(mu)*0.25; % Average chemical potential

chem_elem = mu([hea_vac(:,55) hea_vac(:,1:54)]);

%Generate adjacency matrix
A = zeros(55,55,Nvac); 

for i=1:Nvac
        
    ivac = vacIDs(i);
      
    p0 = ratoms(ivac,2:4);
    
    for l=1:NN4+1
        A(l, l, i) = IntValues(hea_vac(i,l)); %chem_elem(i,l)*(rcut/3)^2;%(chem_elem(i,l))^(2.4)*0.5;
    end

    % Compute diff with the vacancy site
    for k=1:NN4
        kId = neighbours(ivac,k);
        p1 =  ratoms(kId,2:4);
        
        [rij, dx,dy,dz] = get_distance_pbc(p1, p0, Lx, Ly, Lz);
        
         if rij < rcut
             if (rij-b) < tol
                 rstar = rij/b;
             elseif (rij-sqrt(2)*b) < tol
                 rstar = rij/(sqrt(2)*b);
             elseif (rij-sqrt(3)*b) < tol
                 rstar = rij/(sqrt(3)*b);
             elseif (rij-2*b) < tol
                 rstar = rij/(2*b);
             elseif (rij-sqrt(7)*b) < tol
                 rstar = rij/(sqrt(7)*b);
             elseif (rij-3*b) < tol
                 rstar = rij/(3*b);
             elseif (rij-sqrt(12)*b) < tol
                 rstar = rij/(sqrt(12)*b);
             else
                 rstar = rij/(4*b);
             end
           rstar = rij/rcut;

           A(1, k+1, i) = rstar^-1*IntValues(hea_vac(i,k));%rij^2*(chem_elem(i,1)+chem_elem(i,k))*0.5;
           A(k+1, 1, i) = rstar^-1*IntValues(hea_vac(i,k));%rij^2*(chem_elem(i,1)+chem_elem(i,k))*0.5;
         end
    end
    
    % Now check the other atoms
    
    for m=1:NN4
        mId = neighbours(ivac,m);
        p0 = ratoms(mId,2:4);
        
        for n=1:NN4
            nId = neighbours(ivac,n);
            
            if mId ~= nId
                p1 = ratoms(nId,2:4);
                
                [rmn, dx,dy,dz] = get_distance_pbc(p1,p0,Lx,Ly,Lz);
                
                
                if rmn < rcut
                    if (rij-b) < tol
                        rstar = rij/b;
                    elseif (rij-sqrt(2)*b) < tol
                        rstar = rij/(sqrt(2)*b);
                    elseif (rij-sqrt(3)*b) < tol
                        rstar = rij/(sqrt(3)*b);
                    elseif (rij-2*b) < tol
                        rstar = rij/(2*b);
                    elseif (rij-sqrt(7)*b) < tol
                        rstar = rij/(sqrt(7)*b);
                    elseif (rij-3*b) < tol
                        rstar = rij/(3*b);
                    elseif (rij-sqrt(12)*b) < tol
                        rstar = rij/(sqrt(12)*b);
                    else
                        rstar = rij/(4*b);
                    end

                    A(m+1, n+1, i) = rstar^-1*IntValues(hea_vac(i,m));%rmn^2*(chem_elem(i,m+1)+chem_elem(i,n+1))*0.5;
                    A(n+1, m+1, i) = rstar^-1*IntValues(hea_vac(i,n));%rmn^2*(chem_elem(i,m+1)+chem_elem(i,n+1))*0.5;
                end
                
            end
        end
    end
    A(:,:,i) = A(:,:,i)/max(max(A(:,:,i)));
    
end

adjacency = A;

end

