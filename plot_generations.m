function [arg] = plot_generations(Pars,chem_env, ave_chem, sigma, generations, xcoord, ycoord, zcoord)

arg = 0;

% Select what NNS want to investigate 

NNS = 1; % Neighbor shell 1, 2, 3. 

for j=1:Pars.Ngen
    
    for i=1:Pars.Na
        firstNN(i,:) = chem_env(NNS,:,i,j);
    end
    
    max_chem(j,:) = max(firstNN);
    min_chem(j,:) = min(firstNN);
    
    % This is an arbitrary metrics used to define what sample is best!
    xi(j) = sum((ave_chem(NNS,:,j)-0.25).^2)*sum(sigma(NNS,:,j))*sum((max_chem(j,:)-0.25).^2+(min_chem(j,:)-0.25).^2);
    
end

figure(3)
semilogy(1:Pars.Ngen,xi,'o','LineWidth', 1.5)

hold on
plot([0 1000],[0.000185 0.000185],'--r','LineWidth', 1.0)
plot([0 1000],[0.000925 0.000925],'--r','LineWidth', 1.0)

legend ('\xi', 'FontSize', 24)
legend boxoff

xlabel ('Random generation')
ylabel ('Performance metric')

set(gca, 'FontSize', 18, 'LineWidth', 2)

formatSpec = 'Xi_%dx%dx%d.png';
figname = sprintf(formatSpec,Pars.Nx,Pars.Ny,Pars.Nz);

print('-f3',figname,'-dpng')

figure(30)
formatSpec = 'XiTally_%dx%dx%d.png';
figname = sprintf(formatSpec,Pars.Nx,Pars.Ny,Pars.Nz);

histogram(xi,30)
xlabel ('\xi', 'FontSize', 24)
ylabel ('Tally', 'FontSize', 24)
set(gca, 'FontSize', 18, 'LineWidth', 2)
print('-f30',figname,'-dpng')


fprintf('Performance metric successful \n')

%% 

[val,pos]=min(xi);

j =pos;
  
for i=1:Pars.Na
    firstNN(i,:) = chem_env(NNS,:,i,j);
end

max_chem = max(firstNN);
min_chem = min(firstNN);

func = sum((max_chem-0.25).^2+(min_chem-0.25).^2);

% plot mean concentration and std dev

figure(4)
    
    for i=1:Pars.Ne
        errorbar([ 1 2 3 4]',ave_chem(:,i,j),sigma(:,i,j),'o','Color',Pars.colors(i,:),...
            'MarkerSize',10,'LineWidth',2.0,'MarkerFaceColor',Pars.colors(i,:))
        if i == 1
            hold on
        end
    end
    
    hold off
    xticks([1 2 3 4])
    xticklabels({ '1NN','2NN','3NN','4NN'})
    
%     legend ('Fe','Ni','Cr','Co')
%     legend boxoff
    
    xlabel ('Neighbor shell')
    ylabel ('Atomic molar fraction')
    
    grid on
    axis([0.5 4.5 0.05 0.45])
    
    set(gca, 'FontSize', 18, 'LineWidth', 2)
    
print('-f4','FigConcNeighShell.png','-dpng')

%
figure(5)

[hBox hErr] = barwitherr(sigma(:,:,j),ave_chem(:,:,j),'LineWidth',1.50);% Plot with errorbars
axis([0.5 4.5 0.05 0.45])

xticks([1 2 3 4])
xticklabels({ '1NN','2NN','3NN','4NN'})

xlabel('Neighbor shell')
ylabel('Atomic molar fraction')


for i=1:Pars.Ne
    set(hBox(i), 'FaceColor',Pars.colors(i,:),'LineWidth',1.50)
end

set(hErr(:), 'LineWidth', 1.5)

set(gca, 'FontSize', 18, 'LineWidth', 2)

formatSpec = 'ConcNeighShellR%d_%dx%dx%d.png';
figname = sprintf(formatSpec,j,Pars.Nx,Pars.Ny,Pars.Nz);

print('-f5',figname,'-dpng')


%% For validation purposes, plot the configuration and the dist of elements.
conf = randi(Pars.Ngen);

figure(10)
scatter3(xcoord,ycoord,zcoord,100,generations(conf,:),'filled')

figure(20)
histogram(generations(conf,:))

arg = pos;

fprintf('Plot chemical environment successful \n')
end

