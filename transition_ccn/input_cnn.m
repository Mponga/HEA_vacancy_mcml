function [Pars] = input_cnn()

Pars.Nx = 50;                     % Select number of lattice size
Pars.Ny = 50;                     % Select number of lattice size
Pars.Nz = 50;                     % Select number of lattice size
Pars.a0 = 1;                     % Select lattice constant
Pars.latconst = 3.558;                     % Select lattice constant
Pars.b = 0.5*sqrt(2)*Pars.a0;    %Burgers vector
Pars.R0 = 2.1*Pars.b;                 %Cutoff distance for neighbor list
Pars.Ngen = 10;                %Number of random generations
Pars.Ne = 4;                     %Max number of elements
Pars.Na = 0;                     %Number of atoms, will be init in fcc_cluster
Pars.colors = callmycolors();    %Call color pallet for plotting
Pars.Temp = 1500;    %Call color pallet for plotting
Pars.kbT = Pars.Temp*8.617333262145e-5; % Thermodynamic factor in eV
Pars.NN = 12;                     %Nearest neighbor in fcc lattice.

%% How many vacancies
Pars.vac = 10;                   % Predefine the number of vacancies in the sample
Pars.vacid = 0;
Pars.vac_type = zeros(Pars.vac,1);
Pars.actvac = ones(Pars.vac,1);  %Define array to make active vacancy. 

%% Parameters for the convolutional neural network
% 
Pars.pixel_Nx=24;
Pars.pixel_Ny=24;

fprintf('Generation successful : N elements %d, Lx %d, Ly %d, Lz %d\n', Pars.Ne, Pars.Nx, Pars.Ny, Pars.Nz)
fprintf('Vacancy generation successful : Nvac %d\n', Pars.vac)

Pars.max_vac=2.4858; %Max vfe used in the network 
Pars.min_vac=0.8159; %Min vfe used in the network 
Pars.Dvac = Pars.max_vac-Pars.min_vac; 

Pars.maxef = 1.8048127;
Pars.mineb = 6.404611099999999e-01;
Pars.Det   = Pars.maxef - Pars.mineb;


end

