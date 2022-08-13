%% Used to read the neb_output.txt file from the cluster
% This file will generate one figure per vacancy site and try to
% incorporate that as the input of the CNN. 

clear all,clc,close all

model = 'image'; %'graph';%

neb_output_1 = dlmread('neb_output_1.txt'); %could replace neb_output
neb_output_2 = dlmread('neb_output.txt'); %could replace neb_output

neb_output = [neb_output_1; neb_output_2];

colors = callmycolors();    %Call color pallet for plotting

vacIDs    = neb_output(:,2); %The vacancy is now migrated, not in the initial position
vacCenter = neb_output(:,1); %The original position of the vacancy
Eforward  = neb_output(:,9);
Ebackward = neb_output(:,10);

maxef = max(Eforward);
maxeb = max(Ebackward);

minef = min(Eforward);
mineb = min(Ebackward);

Trans_Energy = [(Eforward-minef)/(maxef-minef) (Ebackward-mineb)/(maxeb-mineb)];

%% Read atoms file used to compute the nebs. 

% Check that Nx,Ny Nz correspond to the same as the file. 
[Pars] = input_cnn();

% % This file was modified such that the three spaces between elements was
% % replaced by one element, makes it easier to read the file. 
%make neighbours using the reference and use the relaxed positions later
Atoms = dlmread('readhea-genR10_50x50x50_NEB_255101_NN0.lmp'); 

% The lammps file printed with real coordinates, need to rescale. 
xcoord=Atoms(:,3)/Pars.latconst;
ycoord=Atoms(:,4)/Pars.latconst;
zcoord=Atoms(:,5)/Pars.latconst;
typeID=Atoms(:,2);

%% Read relaxed atoms 
ratoms = dlmread('dump.relax.1.234',' ',9,0);

xmin = .9843791387248487e-01;
xmax = 1.7740156208612734e+02;
ymin = 4.9852784999691774e-01; 
ymax = 1.7740147215000314e+02;
zmin = 4.9756685245696847e-01;
zmax = 1.7740243314754323e+02;

Lx = xmax-xmin;
Ly = ymax-ymin;
Lz = zmax-zmin;

rtypeID = ratoms(:,1);
ratoms(:,2) = ratoms(:,2)*(xmax-xmin) + xmin;
ratoms(:,3) = ratoms(:,3)*(ymax-ymin) + ymin;
ratoms(:,4) = ratoms(:,4)*(zmax-zmin) + zmin;

Pars.Na = length(Atoms);

[Cells] = input_cells([xcoord,ycoord,zcoord], Pars);

[neighbours, NumNeigh] = create_neighbor_list_opt([xcoord ycoord zcoord], Pars, Cells); 

%% 

vac1 = [typeID(neighbours(vacIDs,:))    typeID(vacIDs)];
vac2 = [typeID(neighbours(vacCenter,:)) typeID(vacCenter)];

% vac1 = [flip(typeID(neighbours(vacIDs,1:12))) flip(typeID(neighbours(vacIDs,13:18))) flip(typeID(neighbours(vacIDs,19:42))) flip(typeID(neighbours(vacIDs,43:54))) typeID(vacIDs)];
% imgs = generate_imgs(vac0,Pars);

switch model
    case 'image'
        imgs = generate_double_imgs(vac1/4,vac2/4,Pars);
        Pars.pixel_Nx = 24;
        Pars.pixel_Ny = 24;

    case 'graph'
        adjacency1 = generate_graphs(vac1, Pars, vacIDs, neighbours, typeID, ratoms, Lx, Ly, Lz);
        adjacency2 = generate_graphs(vac2, Pars, vacIDs, neighbours, typeID, ratoms, Lx, Ly, Lz);

        imgs(:,:,1,:) = adjacency1;
        imgs(:,:,2,:) = adjacency2;
        Pars.pixel_Nx = 55;
        Pars.pixel_Ny = 55;

        graphplot =false; 

        if graphplot
            % Plot graph plots for some random vacancies
            rgraph = 6137; %randperm(length(adjacency),1);

            G = graph(adjacency1(:,:,rgraph));

            xgraph = ratoms([vacIDs(rgraph) neighbours(vacIDs(rgraph),:)],2);
            ygraph = ratoms([vacIDs(rgraph) neighbours(vacIDs(rgraph),:)],3);
            zgraph = ratoms([vacIDs(rgraph) neighbours(vacIDs(rgraph),:)],4);

            LWidths = 5*G.Edges.Weight/max(G.Edges.Weight);

            p = plot(G,'XData',xgraph,'YData',ygraph,'ZData',zgraph,'EdgeLabel',G.Edges.Weight,'LineWidth',LWidths)

            p.Marker = 's';
            p.NodeColor = 'r';
            p.MarkerSize = 7;

            % xlabel("Predicted Value")
            % ylabel("True Value")
            title("Adjacency Graph")

            grid on
            % axis([0 1 0 1])

            set(gca, 'FontSize', 18, 'LineWidth', 2)
        end
end


% Plot random figures to check consistency

        idx = randperm(length(imgs),4);

        for i = 1:4
            subplot(1,4,i)
            imshow( imgs(:,:,2,idx(i)) ,'InitialMagnification','fit')
        end
%% Generate the CNN for prediction of the migration energies

% Convolution layer parameters
fs = 2; %Filter size: It seems that the optimal value is 3
nf = 16;%128; %Number of filters: Does not change much beyond 16 or more. 
ps = 2; %Pool size for average region.
do = 0.1; %Dropout probability;
maxepochs = 80; 
miniBatchSize  = 64;
solverName = 'adam';%'rmsprop'; %'sgdm';%'adam'; % GradientDecayFactor
InitialLearnRate = 1e-2; %0.5e-2; %5e-3; %1e-3; %It seems that for larger values the final convergence is better

%% Make Training and Validation sets
%

numObservations = length(vacIDs);

P = [0.70 0.15 0.15];
%P = [0.50 0.25 0.25];

rng(124676512);

idx = randperm(numObservations);

idxTrain = idx(1:round(P(1)*numObservations));
idxValidation = idx(round(P(1)*numObservations)+1:round((P(1)+P(2))*numObservations));
idxTest = idx(round((P(1)+P(2))*numObservations)+1:numObservations);

XTrain = imgs(:,:,:,idxTrain);
YTrain = Trans_Energy(idxTrain,:); 
CTrain = typeID(vacIDs(idxTrain));

XValidation = imgs(:,:,:,idxValidation);
YValidation = Trans_Energy(idxValidation,:);
CValidation = typeID(vacIDs(idxValidation));

XTest = imgs(:,:,:,idxTest);
YTest = Trans_Energy(idxTest,:); 
CTest = typeID(vacIDs(idxTest));

%%
% save('HEA_CNN_input_data.mat','XTrain','YTrain','CTrain','XValidation','YValidation','CValidation','XTest','YTest','CTest',"idx");

%%

%
% Create Network Layers

layers = [
     imageInputLayer([Pars.pixel_Nx Pars.pixel_Ny 2])
    convolution2dLayer(fs,nf,'Padding','same')
    batchNormalizationLayer
    reluLayer
    
    %
    averagePooling2dLayer(ps,'Stride',ps)
    convolution2dLayer(fs,nf,'Padding','same')
    batchNormalizationLayer
    reluLayer
    %

    %
    averagePooling2dLayer(ps,'Stride',ps)
    convolution2dLayer(fs,nf,'Padding','same')
    batchNormalizationLayer
    reluLayer
    %
    
    averagePooling2dLayer(ps,'Stride',ps)
    convolution2dLayer(fs,2*nf,'Padding','same')
    batchNormalizationLayer
    reluLayer
    
    averagePooling2dLayer(ps,'Stride',ps)
    convolution2dLayer(fs,4*nf,'Padding','same')
    batchNormalizationLayer
    reluLayer
    
    convolution2dLayer(fs,4*nf,'Padding','same')
    batchNormalizationLayer
    reluLayer
    
    dropoutLayer(do)
    fullyConnectedLayer(2,'WeightLearnRateFactor',2,'BiasLearnRateFactor',2,'WeightL2Factor',2)
    regressionLayer];

% Train Network

% Minimization solver total iterations is given by:
% validationFrequency * MaxEpochs. 
% To increase # of iterations reduce miniBatchSize or increase MaxEpochs

validationFrequency = floor(numel(YTrain)/miniBatchSize);
options = trainingOptions(solverName, ...
    'MiniBatchSize',miniBatchSize, ...
    'SquaredGradientDecayFactor',0.99,...
    'MaxEpochs',maxepochs, ...
    'InitialLearnRate',InitialLearnRate, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropFactor',0.1, ...
    'LearnRateDropPeriod',20, ...
    'Shuffle','every-epoch', ...
    'ValidationData',{XValidation,YValidation}, ...
    'ValidationFrequency',validationFrequency, ...
    'Plots','training-progress', ...
    'OutputNetwork','best-validation-loss');
%     'Verbose',false);
%     'OutputNetwork','best-validation-loss',...
    %'ValidationData',{XValidation,YValidation}, ...
%     'ValidationData',{GValidation,YValidation}, ...

%%
net = trainNetwork(XTrain,YTrain,layers,options);

trans_net1 = net;
save('Transition_CNN_Network','trans_net1');

% Test Network

YPredicted = predict(net,XValidation);

YTestPredicted = predict(net,XTest);

predictionError = YValidation - YPredicted;

thr = 0.10;
numCorrect = sum(abs(predictionError) < thr);
numValidationImages = numel(YValidation);

accuracy = numCorrect/numValidationImages

% Use the root-mean-square error (RMSE) to measure the differences between the predicted and actual angles of rotation.

squares = predictionError.^2;
rmse = sqrt(mean(squares))

% Compute R2 coefficient
Rsq = 1 - sum((YValidation - YPredicted).^2)/sum((YValidation - mean(YValidation)).^2);

%%
Rsq1 = 1 - sum((YValidation(:,1) - YPredicted(:,1)).^2)/sum((YValidation(:,1) - mean(YValidation(:,1))).^2)
Rsq2 = 1 - sum((YValidation(:,2) - YPredicted(:,2)).^2)/sum((YValidation(:,2) - mean(YValidation(:,2))).^2)
%

figure(10)
dx = 1.5*rmse(1);
plot([0 10], [0 10],'r--')
hold on
plot([0 10], [-dx 10-dx],'r--')
plot([0 10], [ dx 10+dx],'r--')
scatter((YPredicted(:,1)'),(YValidation(:,1)'),50,CValidation,'filled')
xlabel("Predicted Value")
ylabel("True Value")

grid on
axis([0 1 0 1])

set(gca, 'FontSize', 18, 'LineWidth', 2)

text(.1,.75,'Forward','FontSize', 18)
txt = ['RMSE ',num2str(rmse(1),'%0.2f')];
txt2 = ['R^2',num2str(Rsq1,'%0.2f')];
text(.1,.7,txt,'FontSize', 18)
text(.1,.65,txt2,'FontSize', 18)
print('-f10','FigForward.png','-dpng')

figure(20)
dx = 1.5*rmse(2);
plot([0 10], [0 10],'r--')
hold on
plot([0 10], [-dx 10-dx],'r--')
plot([0 10], [ dx 10+dx],'r--')
scatter((YPredicted(:,2)'),(YValidation(:,2)'),50,CValidation,'filled')

xlabel("Predicted Value")
ylabel("True Value")

grid on
axis([0 1 0 1])

set(gca, 'FontSize', 18, 'LineWidth', 2)

text(.1,.75,'Backward','FontSize', 18)
txt = ['RMSE ',num2str(rmse(2),'%0.2f')];
txt2 = ['R^2',num2str(Rsq2,'%0.2f')];

text(.1,.7,txt,'FontSize', 18)
text(.1,.65,txt2,'FontSize', 18)

print('-f20','FigBackward.png','-dpng')

figure(30)
dx = 1.5*rmse(2);
plot([0 10], [0 10],'r--')
hold on
plot([0 10], [-dx 10-dx],'r--')
plot([0 10], [ dx 10+dx],'r--')
scatter((YTestPredicted(:,1)'),(YTest(:,1)'),50,CTest,'filled')
hold on
scatter((YTestPredicted(:,2)'),(YTest(:,2)'),50,CTest,'filled')
xlabel("Predicted Value")
ylabel("True Value")
title("Test Set")

grid on
axis([0 1 0 1])

set(gca, 'FontSize', 18, 'LineWidth', 2)

text(.1,.75,'Forward and Backward','FontSize', 18)
txt = ['RMSE ',num2str(rmse(2),'%0.2f')];
txt2 = ['R^2',num2str(Rsq1,'%0.2f')];

text(.1,.7,txt,'FontSize', 18)
text(.1,.65,txt2,'FontSize', 18)

print('-f30','FigTestSet.png','-dpng')
