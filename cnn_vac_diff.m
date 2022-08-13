%% This program generates a fcc random solid solution hea with Ne
%  It also generates the neighbor list for each atom (in parallel) and
%  computes the chemical environment of the sample. In addition to that, it
%  also computes the images for convolutional neural network evaluation of
%  the vacancy formation and migration energy. All together, it also has
%  the ability to perform Monte Carlo simulation of vacancy diffusion to investigate the paths for diffusion. 

clear all, clc, close all 

sample = 'new'; %'large'
chemenv = 'no'; %'yes'

colors = callmycolors();    %Call color pallet for plotting

[Pars] = input_cnn();

[xcoord,ycoord,zcoord,Pars.Na]=fcc_cluster(Pars);

[Cells] = input_cells([xcoord,ycoord,zcoord], Pars);

% Make several mixtures to assess the chem environment
generations = make_random_generations(Pars);

Pars.vac_id = randsample(Pars.Na,Pars.vac);

Nflips=1e6;
% Create neighbor list for all atoms.

[neighbours, NumNeigh] = create_neighbor_list_opt([xcoord ycoord zcoord], Pars, Cells); 

% Analyze the chemical environment of the sample 
switch chemenv
    case 'yes'
        % Plot the chemical environment
        [chem_env, ave_chem, sigma] = analyze_chem_env(Pars,neighbours,NumNeigh, generations);
        
        % plot the one with minimum R^2 to equiatomic in the 1NN
        best_gen=plot_generations(Pars,chem_env, ave_chem, sigma, generations, xcoord, ycoord, zcoord);
        
    case 'no'
        best_gen = 1;
end

% Make a 4D array for each image [Nx_pixels x Ny_pixels x 1 x Pars.Na]
%  Generate a 4D-double array with the figures for each site.
mygen = generations(best_gen,:); %id of the random generated sites.

Pars.vac_type = mygen(Pars.vac_id)';

fprintf('%d Vacancies with atom ID and type\n',Pars.vac);

for i=1:Pars.vac
    fprintf('Vac #%d: %d %d \n',i, Pars.vac_id(i), Pars.vac_type(i));
end
fprintf('\n');

hea_ids = zeros(Pars.Na,NumNeigh(1)+1);

for i=1:Pars.Na
    hea_ids(i,:) = [mygen(neighbours(i,:)) mygen(i)];
end

imgs = generate_large_imgs(hea_ids,Pars); %generate_imgs(hea_ids,Pars);

% Plot random figures to check consistency

idx = randperm(Pars.Na,4);

for i = 1:numel(idx)
    subplot(1,4,i)       
    imshow( imgs(:,:,:,idx(i)) ,'InitialMagnification','fit')
end

% Laod the previously save network
load Vac_CNN_Network.mat

% YPredicted = (Pars.Dvac)*predict(vacnet1,imgs(:,:,:,1:100:end))+Pars.min_vac; %rescale values. 
% histogram(YPredicted);
% vfe_mean = mean(YPredicted);

% Now load Transition network
load Transition_CNN_Network.mat

%% One could restart the problem by selecting new vacancies. 
clc, close all

mytemp = [300:150:1500];

for i=1:length(mytemp)

    mywaitbar = 'no'; %'yes'

    Pars.Temp = mytemp(i); % Thermodynamic factor in eV
    Pars.kbT = Pars. Temp*8.617333262145e-5; % Thermodyn0amic factor in eV

    cases = 500;
    limit0 = 1e-6; %time limit for the diffusion simulation. (1 us)

    limit = limit0*exp(-(Pars.Temp-300)/(300));

    Dv = zeros(cases,1);

    % Need to include some sort of critical rate for time update since there
    % are some very fast processes that get the diffusion stack and time does
    % not move. These are fa st processes involving quick swap of vacancies.
    if Pars.Temp == 300
        rate_critical = 700;
    elseif Pars.Temp == 450
        rate_critical = 8500;
    elseif Pars.Temp == 600
        rate_critical = Nflips/(2*limit*1e9);%1000;
    elseif Pars.Temp == 900
        rate_critical = 5000;
    elseif Pars.Temp == 1200
        rate_critical = 15000;
    elseif Pars.Temp == 1500
        rate_critical = 25000;
    else
        rate_critical = 25000;
    end

    filename = sprintf('diffusion_%dK.txt', Pars.Temp );
    fileID = fopen(filename,'a');

     fprintf(fileID,'Starting diffusion simulations\n');
     fprintf(fileID,'Number of Cases: %d\n',cases);
     fprintf(fileID,'Temperature: %d\n',Pars.Temp);
     fprintf(fileID,'Maximum Flips: %d\n',Nflips);
     fprintf(fileID,'Limit Time: %12.1f [ns] \n',limit*1e9);
     fprintf(fileID,'Critical rate: %12.1f [it/ns] \n',rate_critical);

    fprintf('Starting diffusion simulations\n');
    fprintf('Number of Cases: %d\n',cases);
    fprintf('Temperature: %d\n',Pars.Temp);
    fprintf('Maximum Flips: %d\n',Nflips);
    fprintf('Limit Time: %12.1f [ns] \n',limit*1e9);
    fprintf('Critical rate: %12.1f [it/ns] \n',rate_critical);

    for l = 1:cases %number of random generated vacancy cases to benchmark diffusive behavior

        Pars.vac = 1;                   % Predefine the number of vacancies in the sample
        Pars.vac_id = randsample(Pars.Na,Pars.vac);
        Pars.vac_type = mygen(Pars.vac_id)';

        fprintf(fileID,'%d %d %d ', l, Pars.vac_id, Pars.vac_type);
        fprintf('%d %d %d ',l, Pars.vac_id, Pars.vac_type)

        % Loop Nflips times to move vacancies

        vac_IDs = zeros(Pars.vac,Nflips+1);
        vac_IDs(:,1) = Pars.vac_id;

        vfe_val = zeros(Pars.vac,Nflips+1);

        vac0 = [mygen(neighbours(Pars.vac_id,:)) Pars.vac_type];

        new_img = generate_large_imgs(vac0,Pars); %Generate nvac new images only.

        vfe_val(:,1) = (Pars.Dvac)*predict(vacnet1,new_img)+Pars.min_vac; %rescale values.

        time = zeros(Pars.vac);
        mytime = zeros(Nflips,1);
        vac_distance = zeros(Nflips,1);

        Nflips = 20000;
        switch mywaitbar
            case 'yes'
                % Create waitbar
                fwaitbar = waitbar(0,'1','Name','Computing vacancy diffusion',...
                    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

                setappdata(fwaitbar,'canceling',0);
        end
        %

        for i=1:Nflips

            switch mywaitbar
                case 'yes'
                    if getappdata(fwaitbar,'canceling')
                        break
                        delete(fwaitbar)
                    end
            end

            Pars.actvac = check_vac([xcoord ycoord zcoord], Pars); % Check if vac are close enough or not.

            % Look over all NN1 and compute energy of new sites, then select neighbor to flip with the vacancy
            [flipval, swap, kt] = swap_vacancies(Pars,vacnet1,trans_net1,mygen,neighbours,Pars.Dvac,Pars.min_vac,vfe_val(:,i));

            delta_t  = log(1/rand(1))*(1/(Pars.nu0*kt)); %-log(rand(1))/kt;
            time = time + delta_t;
            mytime(i) = time*1e9; %time in nanoseconds

            %     fprintf('%8.4e\n', time);

            % Flip only active vacancies
            actswaps = Pars.actvac.*swap;
            myswaps = nonzeros(actswaps);
            positions = length(myswaps)/2;
            swap(1:positions,1) = myswaps(1:positions);
            swap(1:positions,2) = myswaps(positions+1:length(myswaps));

            % Once the best sites are identified, make the swap between site and vacancy
            Pars.vac_id = swap(:,1); %Change vac ID only
            mygen(swap) = mygen(flip(swap,2)); %c hange atoms type as mesh does not move

            vac_IDs(:,i+1) = Pars.vac_id;
            vfe_val(:,i+1) = flipval;

            if i==1
                myvacs = vac_IDs(1);
                p0 = [xcoord(myvacs) ycoord(myvacs) zcoord(myvacs)];
            end
            myvacs = vac_IDs(:,i);
            p1 = [xcoord(myvacs) ycoord(myvacs) zcoord(myvacs)];

            [R, dx, dy, dz] = get_distance_pbc(p1,p0,Pars.Nx,Pars.Ny,Pars.Nz);

            vac_distance(i) = R*Pars.b*Pars.latconst*1e-10; %distance in m

            if time > limit && i > 100
                break
            end

            rate = i/(time*1e9); %iterations/ns

            if rem(i,200)==0
                if rate > rate_critical
                    %l = l-1;
                    fprintf(fileID,'failed \n');
                    fprintf('failed \n');
                    break
                end
            end
            switch mywaitbar
                case 'yes'
                    % Update waitbar and message
                    waitbar(time/limit,fwaitbar,sprintf('Case: %d - %d/%d %12.1f it/ns %12.1f [ns] out of %12.1f [ns]',l, i, Nflips, rate, time*1e9, limit*1e9))
                    %
            end

        end

        switch mywaitbar
            case 'yes'
                delete(fwaitbar)
        end

        if i == Nflips
            l = l-1;
            fprintf('failed \n');
            fprintf(fileID,'failed \n');
        else
            vac_movmean = movmean(vac_distance,14);

            % plot(vac_movmean,'LineWidth',2)
            %         semilogx(mytime/1e9,vac_movmean,'LineWidth',2)
            semilogx(mytime/1e9,vac_distance*1e9,'LineWidth',2)
            hold on
            ylabel('<(r-r(0))^2> [nm]')
            xlabel('Time [s]')
            set(gca, 'FontSize', 18, 'LineWidth', 2)

            drawnow


            r_ave_2 = mean(nonzeros(vac_distance).^2);

            Dv(l) = r_ave_2/(6*time);

            fprintf(fileID,'%8.4e %8.4e %8.4e\n', sqrt(r_ave_2), time, Dv(l));
            fprintf('%8.4e %8.4e %8.4e\n', sqrt(r_ave_2), time, Dv(l))
        end
    end

    fclose(fileID);

end
%%
myplots = 0; 

if myplots == 1
    close all
    figure('Renderer', 'painters', 'Position', [10 10 900 600])

    subplot(1,2,1)

    for i=1:Pars.vac

        plot(vfe_val(i,:),'Color',Pars.colors(i,:),'LineWidth',1.50)
        hold on
    end

    % Plot spatially the vacancies over time

    subplot(1,2,2)
    myvacs = vac_IDs(:,1);
    scatter3(xcoord(myvacs),ycoord(myvacs),zcoord(myvacs),100,Pars.vac_type,'filled')
    axis([0 Pars.Nx 0 Pars.Ny 0 Pars.Nz])

    for i=1:Nflips

        formatSpec = 'Iteration %d - Ellapsed time %8.2f [ns]';
        text = sprintf(formatSpec,i,mytime(i));
        title(text)
        myvacs = vac_IDs(:,i);
        scatter3(xcoord(myvacs),ycoord(myvacs),zcoord(myvacs),100,Pars.vac_type,'filled')
        hold on
        axis([0 Pars.Nx 0 Pars.Ny 0 Pars.Nz])
        pbaspect([1 1 1])
        drawnow

    end

end

%% Plot results for manuscript

experiments = [1073 7e-19
       1173 9e-18
       1373 2e-16];

data_300K = dlmread('diffusion_300K.txt',' ',6,0);
data_450K = dlmread('diffusion_450K.txt',' ',6,0);
data_600K = dlmread('diffusion_600K.txt',' ',6,0);
data_750K = dlmread('diffusion_750K.txt',' ',6,0);
data_900K = dlmread('diffusion_900K.txt',' ',6,0);
data_1050K = dlmread('diffusion_1050K.txt',' ',6,0);
data_1200K = dlmread('diffusion_1200K.txt',' ',6,0);
data_1350K = dlmread('diffusion_1350K.txt',' ',6,0);
data_1500K = dlmread('diffusion_1500K.txt',' ',6,0);


colors = callmycolors();    %Call color pallet for plotting

Pars.kb = 8.617333262145e-5; % Thermodynamic factor in eV

median_300K = median(data_300K(:,6));
median_450K = median(data_450K(:,6));
median_600K = median(data_600K(:,6));
median_750K = median(data_750K(:,6));
median_900K = median(data_900K(:,6));
median_1050K = median(data_1050K(:,6));
median_1200K = median(data_1200K(:,6));
median_1200K = median(data_1350K(:,6));
median_1500K = median(data_1500K(:,6));

mean_300K = mean(data_300K(:,6));
mean_450K = mean(data_450K(:,6));
mean_600K = mean(data_600K(:,6));
mean_750K = mean(data_750K(:,6));
mean_900K = mean(data_900K(:,6));
mean_1200K = mean(data_1200K(:,6));
mean_1050K = mean(data_1050K(:,6));
mean_1050K = mean(data_1350K(:,6));
mean_1500K = mean(data_1500K(:,6));

close all

h = zeros(1,4);

figure(1)

h(1) = semilogy(1,1,'o','Color',colors(1,:),'MarkerFaceColor',colors(1,:),'DisplayName','Ni')
hold on
h(2) = semilogy(1,1,'s','Color',colors(2,:),'MarkerFaceColor',colors(2,:),'DisplayName','Fe')
h(3) = semilogy(1,1,'^','Color',colors(3,:),'MarkerFaceColor',colors(3,:),'DisplayName','Cr')
h(4) = semilogy(1,1,'v','Color',colors(4,:),'MarkerFaceColor',colors(4,:),'DisplayName','Co')


semilogy(10000/300,data_300K(data_300K(:,3) == 1,6),'o','Color',colors(1,:),'MarkerFaceColor',colors(1,:))
hold on
semilogy(10000/300,data_300K(data_300K(:,3) == 2,6),'s','Color',colors(2,:),'MarkerFaceColor',colors(2,:))
semilogy(10000/300,data_300K(data_300K(:,3) == 3,6),'^','Color',colors(3,:),'MarkerFaceColor',colors(3,:))
semilogy(10000/300,data_300K(data_300K(:,3) == 4,6),'v','Color',colors(4,:),'MarkerFaceColor',colors(4,:))


semilogy(10000/450,data_450K(data_450K(:,3) == 1,6),'o','Color',colors(1,:),'MarkerFaceColor',colors(1,:))
hold on
semilogy(10000/450,data_450K(data_450K(:,3) == 2,6),'s','Color',colors(2,:),'MarkerFaceColor',colors(2,:))
semilogy(10000/450,data_450K(data_450K(:,3) == 3,6),'^','Color',colors(3,:),'MarkerFaceColor',colors(3,:))
semilogy(10000/450,data_450K(data_450K(:,3) == 4,6),'v','Color',colors(4,:),'MarkerFaceColor',colors(4,:))

semilogy(10000/600,data_600K(data_600K(:,3) == 1,6),'o','Color',colors(1,:),'MarkerFaceColor',colors(1,:))
hold on
semilogy(10000/600,data_600K(data_600K(:,3) == 2,6),'s','Color',colors(2,:),'MarkerFaceColor',colors(2,:))
semilogy(10000/600,data_600K(data_600K(:,3) == 3,6),'^','Color',colors(3,:),'MarkerFaceColor',colors(3,:))
semilogy(10000/600,data_600K(data_600K(:,3) == 4,6),'v','Color',colors(4,:),'MarkerFaceColor',colors(4,:))

semilogy(10000/750,data_750K(data_750K(:,3) == 1,6),'o','Color',colors(1,:),'MarkerFaceColor',colors(1,:))
hold on
semilogy(10000/750,data_750K(data_750K(:,3) == 2,6),'s','Color',colors(2,:),'MarkerFaceColor',colors(2,:))
semilogy(10000/750,data_750K(data_750K(:,3) == 3,6),'^','Color',colors(3,:),'MarkerFaceColor',colors(3,:))
semilogy(10000/750,data_750K(data_750K(:,3) == 4,6),'v','Color',colors(4,:),'MarkerFaceColor',colors(4,:))

semilogy(10000/900,data_900K(data_900K(:,3) == 1,6),'o','Color',colors(1,:),'MarkerFaceColor',colors(1,:))
hold on
semilogy(10000/900,data_900K(data_900K(:,3) == 2,6),'s','Color',colors(2,:),'MarkerFaceColor',colors(2,:))
semilogy(10000/900,data_900K(data_900K(:,3) == 3,6),'^','Color',colors(3,:),'MarkerFaceColor',colors(3,:))
semilogy(10000/900,data_900K(data_900K(:,3) == 4,6),'v','Color',colors(4,:),'MarkerFaceColor',colors(4,:))

semilogy(10000/1050,data_1050K(data_1050K(:,3) == 1,6),'o','Color',colors(1,:),'MarkerFaceColor',colors(1,:))
hold on
semilogy(10000/1050,data_1050K(data_1050K(:,3) == 2,6),'s','Color',colors(2,:),'MarkerFaceColor',colors(2,:))
semilogy(10000/1050,data_1050K(data_1050K(:,3) == 3,6),'^','Color',colors(3,:),'MarkerFaceColor',colors(3,:))
semilogy(10000/1050,data_1050K(data_1050K(:,3) == 4,6),'v','Color',colors(4,:),'MarkerFaceColor',colors(4,:))

semilogy(10000/1200,data_1200K(data_1200K(:,3) == 1,6),'o','Color',colors(1,:),'MarkerFaceColor',colors(1,:))
hold on
semilogy(10000/1200,data_1200K(data_1200K(:,3) == 2,6),'s','Color',colors(2,:),'MarkerFaceColor',colors(2,:))
semilogy(10000/1200,data_1200K(data_1200K(:,3) == 3,6),'^','Color',colors(3,:),'MarkerFaceColor',colors(3,:))
semilogy(10000/1200,data_1200K(data_1200K(:,3) == 4,6),'v','Color',colors(4,:),'MarkerFaceColor',colors(4,:))

semilogy(10000/1350,data_1350K(data_1350K(:,3) == 1,6),'o','Color',colors(1,:),'MarkerFaceColor',colors(1,:))
hold on
semilogy(10000/1350,data_1350K(data_1350K(:,3) == 2,6),'s','Color',colors(2,:),'MarkerFaceColor',colors(2,:))
semilogy(10000/1350,data_1350K(data_1350K(:,3) == 3,6),'^','Color',colors(3,:),'MarkerFaceColor',colors(3,:))
semilogy(10000/1350,data_1350K(data_1350K(:,3) == 4,6),'v','Color',colors(4,:),'MarkerFaceColor',colors(4,:))

semilogy(10000/1500,data_1500K(data_1500K(:,3) == 1,6),'o','Color',colors(1,:),'MarkerFaceColor',colors(1,:))
hold on
semilogy(10000/1500,data_1500K(data_1500K(:,3) == 2,6),'s','Color',colors(2,:),'MarkerFaceColor',colors(2,:))
semilogy(10000/1500,data_1500K(data_1500K(:,3) == 3,6),'^','Color',colors(3,:),'MarkerFaceColor',colors(3,:))
semilogy(10000/1500,data_1500K(data_1500K(:,3) == 4,6),'v','Color',colors(4,:),'MarkerFaceColor',colors(4,:))

D01 = 0.1e-5;
D02 = 0.001e-5;

T=300:100:1500;

semilogy(10000./T,D01*exp(-1.9./(Pars.kb.*T)),'k--','LineWidth',2.0)

semilogy(10000./T,D02*exp(-2.7./(Pars.kb.*T)),'k--','LineWidth',2.0)

semilogy(10000./T,0.1e-6*exp(-2.3./(Pars.kb.*T)),'r-','LineWidth',2.0)

semilogy(10000./experiments(:,1),experiments(:,2),'d','Color','b','MarkerFaceColor',colors(6,:),'MarkerSize',10)

axis([5 25 1e-40 1e-10])
legend(h(1:4));
legend boxoff

ylabel('Diffusivity [m^2/s]')
xlabel('10000/T [1/K]')
set(gca, 'FontSize', 18, 'LineWidth', 2.5)
grid on

h=text(17.5,5e-22,'E_f = 1.5 eV','FontSize',14)
set(h,'Rotation',-25);
h=text(17.5,1e-26,'Median','FontSize',14)
set(h,'Rotation',-25);
h=text(17.5,5e-31,'E_f = 2.7 eV','FontSize',14)
set(h,'Rotation',-35);

formatSpec = 'Diffusivity.png';
figname = sprintf(formatSpec);
print('-f1',figname,'-dpng')
