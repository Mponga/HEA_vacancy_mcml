%% Make a training and validation set to train a regresion model using CNN
% Include the distance in this version to try to make it more accurate. 
% This script seems to be the most accurate net for predicting vacancies and interstitials.
% Revised on June 3, 2022.

clear all,clc, close all

defect =  'vacancy'; %'interstitial';
element = 'All'; %'Ni'; %Cr Fe Ni
mydata = 'simple';%'simple'; %augmented data for better recognition pattern. 

Na = 4; % # of atoms type in sample
vfe_col =56; % column where vfe is located

% Input layer size
Nx = 24;%12; % Pixels in x-dirextion
Ny = 24;%12; % Pixels in y-dirextion

% Convolution layer parameters
fs = 2;%3; %Filter size: It seems that the optimal value is 3
nf = 16;%32; %Number of filters: Does not change much beyond 16 or more. 
ps = 2; %Pool size for average region.
do = 0.1; %Dropout probability;
maxepochs = 45; 
miniBatchSize  = 32;    

switch element
    
    case 'Co'
        
        switch defect
            case 'vacancy'
                hea_vac = dlmread('NN4_Co_Vfe.txt',' ',2,0);
            case 'interstitial'
                hea_vac = dlmread('NN4_Co_Sfe.txt',' ',2,0);
        end   
        
    case 'Cr'
        
        switch defect
            case 'vacancy'
                hea_vac = dlmread('NN4_Cr_Vfe.txt',' ',2,0);
            case 'interstitial'
                hea_vac = dlmread('NN4_Cr_Sfe.txt',' ',2,0);
        end 
        
    case 'Fe'
        
        switch defect
            case 'vacancy'
                hea_vac = dlmread('NN4_Fe_Vfe.txt',' ',2,0);
            case 'interstitial'
                hea_vac = dlmread('NN4_Fe_Sfe.txt',' ',2,0);
        end 
        
    case 'Ni'
        
        switch defect
            case 'vacancy'
                hea_vac = dlmread('NN4_Ni_Vfe.txt',' ',2,0);
            case 'interstitial'
                hea_vac = dlmread('NN4_Ni_Sfe.txt',' ',2,0);
        end 
        
    case 'All'
                   
        switch defect
            case 'vacancy'
                hea_vac_Co = dlmread('NN4_Co_Vfe.txt',' ',2,0);
                hea_vac_Cr = dlmread('NN4_Cr_Vfe.txt',' ',2,0);
                hea_vac_Fe = dlmread('NN4_Fe_Vfe.txt',' ',2,0);
                hea_vac_Ni = dlmread('NN4_Ni_Vfe.txt',' ',2,0);
            case 'interstitial'
                hea_vac_Co = dlmread('NN4_Co_Sfe.txt',' ',2,0);
                hea_vac_Cr = dlmread('NN4_Cr_Sfe.txt',' ',2,0);
                hea_vac_Fe = dlmread('NN4_Fe_Sfe.txt',' ',2,0);
                hea_vac_Ni = dlmread('NN4_Ni_Sfe.txt',' ',2,0);
        end
       
        max_vac = max([ max(hea_vac_Co(:,vfe_col)) max(hea_vac_Cr(:,vfe_col)) max(hea_vac_Fe(:,vfe_col)) max(hea_vac_Ni(:,vfe_col))]);
        min_vac = min([ min(hea_vac_Co(:,vfe_col)) min(hea_vac_Cr(:,vfe_col)) min(hea_vac_Fe(:,vfe_col)) min(hea_vac_Ni(:,vfe_col))]);
                
        hea_vac_Co(:,vfe_col) = (hea_vac_Co(:,vfe_col)-min_vac)/(max_vac-min_vac); %Make 0->1
        hea_vac_Cr(:,vfe_col) = (hea_vac_Cr(:,vfe_col)-min_vac)/(max_vac-min_vac); %Make 0->1
        hea_vac_Fe(:,vfe_col) = (hea_vac_Fe(:,vfe_col)-min_vac)/(max_vac-min_vac); %Make 0->1
        hea_vac_Ni(:,vfe_col) = (hea_vac_Ni(:,vfe_col)-min_vac)/(max_vac-min_vac); %Make 0->1

        hea_vac = [hea_vac_Co; hea_vac_Cr; hea_vac_Fe; hea_vac_Ni];
        
end

%% Make augmented data using only the train set
%  Augmented data 
switch mydata
    case 'augmented'
        NN1 = 12;
        NN2 = 18;
        NN3 = 42;
        NN4 = 54;
        
         reverse_idx = [flip(1:NN1) flip(NN1+1:NN2) flip(NN2+1:NN3) flip(NN3+1:NN4) 55 56];
%         
%          reverse2_idx = [flip(1:6) flip(7:12) flip(13:15)  flip(16:18) flip(19:30) flip(31:42) flip(43:48) flip(49:54) 55 56];
%                 
%         new_hea_vac = [hea_vac(:,reverse_idx); hea_vac(:,reverse2_idx)];
%         
%         hea_vac = [hea_vac; new_hea_vac];
        
%         reverse_idx = [flip(1:6) flip(7:12) flip(13:15)  flip(16:18) flip(19:30) flip(31:42) flip(43:48) flip(49:54) 55 56];

        new_hea_vac = hea_vac(:,reverse_idx); 
        hea_vac = [hea_vac; new_hea_vac];

%         hea_vac = new_hea_vac;
     
                
        fprintf('\n\n Using data augmentation - %d images\n\n',length(hea_vac));

    case 'simple'
        fprintf('\n\n No data augmentation - Using %d images\n\n',length(hea_vac));
end

%%

Nsamples  = length(hea_vac);
hea_vac = hea_vac(1:Nsamples,:);

%% Dimensionalize the values such that all components go 0->1
% 
hea_vac(:,1:vfe_col-1) = hea_vac(:,1:vfe_col-1)/Na; %Make gray scale
% 
NN1_conc = [hea_vac(:,1:12)] *4;

% Counts how many elements are in the first neighbor shell
for i=1:length(NN1_conc)
    
    A = NN1_conc(i,:);
    c = unique(A);
    
    for j = 1:length(c)
        counts(i,j) = sum(A==c(j)); % number of times each unique value is repeated
    end
    
    [myval,myidx]=max(counts(i,:));
    
    Val(i,1) = myval;
    Idx(i,1) = myidx;
    
end

%% Make a 4D array for each image [4 x 5 x 1 x samples in hea_vac]
%  Generate a 4D-double array with the figures and color layer
% New image  format
newimgs = zeros(Nx,Ny,1,length(hea_vac));

imgs = generate_large_imgs(hea_vac,Nx,Ny);
%% Make Training and Validation sets

numObservations = length(hea_vac);
P = [0.70 0.15 0.15];
rng(124676512);
idx = randperm(numObservations);

idxTrain = idx(1:round(P(1)*numObservations));
idxValidation = idx(round(P(1)*numObservations)+1:round((P(1)+P(2))*numObservations));
idxTest = idx(round((P(1)+P(2))*numObservations)+1:numObservations);

XTrain = imgs(:,:,:,idxTrain);
YTrain = hea_vac(idxTrain,vfe_col); 
CTrain = Idx(idxTrain);

XValidation = imgs(:,:,:,idxValidation);
YValidation = hea_vac(idxValidation,vfe_col);
CValidation = Idx(idxValidation);
VValidation = Val(idxValidation,:);

XTest = imgs(:,:,:,idxTest);
YTest = hea_vac(idxTest,vfe_col); 
CTest = Idx(idxTest);

%% Plot random figures to check consistency

numTrainImages = numel(YTrain);
figure

idx = randperm(numTrainImages,4);

for i = 1:numel(idx)
    subplot(1,4,i)       
    imshow( XTrain(:,:,:,idx(i)) ,'InitialMagnification','fit')
end

%% Create Network Layers

layers = [
    imageInputLayer([Nx Ny 1])
    convolution2dLayer(fs,nf,'Padding','same')
    batchNormalizationLayer
    reluLayer

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
    fullyConnectedLayer(1)
    regressionLayer];

%% Train Network

% Minimization solver total iterations is given by:
% validationFrequency * MaxEpochs. 
% To increase # of iterations reduce miniBatchSize or increase MaxEpochs

solverName = 'adam';%'rmsprop'; %'sgdm';%'adam';
% GradientDecayFactor
InitialLearnRate = 0.5e-2;% 5e-2;  0.5e-2; %1e-3; %It seems that for larger values the final convergence is better

validationFrequency = floor(numel(YTrain)/miniBatchSize);
options = trainingOptions(solverName, ...
    'SquaredGradientDecayFactor',0.99,...
    'MiniBatchSize',miniBatchSize, ...
    'MaxEpochs',maxepochs, ...
    'InitialLearnRate',InitialLearnRate, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropFactor',0.1, ...
    'LearnRateDropPeriod',20, ...
    'Shuffle','every-epoch', ...
    'ValidationData',{XValidation,YValidation}, ...
    'ValidationFrequency',validationFrequency, ...
    'Plots','training-progress', ...
    'OutputNetwork','best-validation-loss');%, ...
    %'Verbose',false);

net = trainNetwork(XTrain,YTrain,layers,options);

vacnet1 = net;
save('Vac_CNN_Network','vacnet1');

%% Test Network

YPredicted = predict(net,XValidation);

predictionError = YValidation - YPredicted;

thr = 0.05;
numCorrect = sum(abs(predictionError) < thr);
numValidationImages = numel(YValidation);

accuracy = numCorrect/numValidationImages

% Use the root-mean-square error (RMSE) to measure the differences between the predicted and actual angles of rotation.

squares = predictionError.^2;
rmse = sqrt(mean(squares))

% Compute R2 coefficient

Rsq = 1 - sum((YValidation - YPredicted).^2)/sum((YValidation - mean(YValidation)).^2)

%% Visualize Predictions
% Visualize the predictions in a scatter plot. Plot the predicted values against the true values.
%

plot_type = 'dimensionless';%'dimensionless'; %'physical';%'dimensionless'; %'physical';%'dimensionless'; %'physical';

figure(3)

switch plot_type
    case 'dimensionless'
        dx = 1.5*rmse;
        plot([0 10], [0 10],'r--')
        hold on
        plot([0 10], [-dx 10-dx],'r--')
        plot([0 10], [ dx 10+dx],'r--')
        
        scatter((YPredicted'),(YValidation'),10*VValidation',CValidation','filled')
        
        xlabel("Predicted Value")
        ylabel("True Value")
        
        grid on
        axis([0 1 0 1])
        
        set(gca, 'FontSize', 18, 'LineWidth', 2)

        switch element
            case 'Co'
                text(.225,.75,'Co Vacancies','FontSize', 18)
                txt = ['RMSE ',num2str(rmse,'%f')];
                text(.2,.7,txt,'FontSize', 18)
                print('-f3','FigCoRegression.png','-dpng')
            case 'Cr'
                text(.225,.75,'Cr Vacancies','FontSize', 18)
                txt = ['RMSE ',num2str(rmse,'%f')];
                text(.2,.7,txt,'FontSize', 18)
                print('-f3','FigCrRegression.png','-dpng')
            case 'Fe'
                text(.225,.75,'Fe Vacancies','FontSize', 18)
                txt = ['RMSE ',num2str(rmse,'%f')];
                text(.2,.7,txt,'FontSize', 18)
                print('-f3','FigFeRegression.png','-dpng')
            case 'Ni'
                text(.225,.75,'Ni Vacancies','FontSize', 18)
                txt = ['RMSE ',num2str(rmse,'%f')];
                text(.2,.7,txt,'FontSize', 18)
                print('-f3','FigNiRegression.png','-dpng')
            case 'All'
                text(.225,.75,'All SIAs','FontSize', 18)
                txt = ['RMSE ',num2str(rmse,'%.03f') ,' & R^2 ',num2str(Rsq,'%.04f')];
                text(.1,.8,txt,'FontSize', 18)
                print('-f3','FigAllRegression.png','-dpng')
        end
        
    case 'physical'
        
        dx = rmse*(mvf-minvf);
        plot([0 10], [0 10],'r--')
        hold on
        plot([0 10], [-dx 10-dx],'r--')
        plot([0 10], [ dx 10+dx],'r--')
        
        scatter((YPredicted')*(mvf-minvf)+minvf,(YValidation')*(mvf-minvf)+minvf,10*VValidation',CValidation','filled')
        xlabel("Predicted Value")
        ylabel("True Value")
        
        grid on
        axis([minvf mvf minvf mvf])
        set(gca, 'FontSize', 18, 'LineWidth', 2)

        
        switch element
            case 'Co'
                text(.35*mvf,.75*mvf,'Co Vacancies','FontSize', 18)
                txt = ['RMSE ',num2str(rmse,'%f')];
                text(.35*mvf,.7*mvf,txt,'FontSize', 18)
                print('-f3','FigCoRegression.png','-dpng')
            case 'Cr'
                text(.35*mvf,.75*mvf,'Cr Vacancies','FontSize', 18)
                txt = ['RMSE ',num2str(rmse,'%f')];
                text(.35*mvf,.7*mvf,txt,'FontSize', 18)
                print('-f3','FigCrRegression.png','-dpng')
            case 'Fe'
                text(.35*mvf,.75*mvf,'Fe Vacancies','FontSize', 18)
                txt = ['RMSE ',num2str(rmse,'%f')];
                text(.35*mvf,.7*mvf,txt,'FontSize', 18)
                print('-f3','FigFeRegression.png','-dpng')
            case 'Ni'
                text(.35*mvf,.75*mvf,'Ni Vacancies','FontSize', 18)
                txt = ['RMSE ',num2str(rmse,'%f')];
                text(.35*mvf,.7*mvf,txt,'FontSize', 18)
                print('-f3','FigNiRegression.png','-dpng')
            case 'All'
                text(.35*mvf,.75*mvf,'All SIAs','FontSize', 18)
                txt = ['RMSE ',num2str(rmse,'%.03f') ,' & R^2 ',num2str(Rsq,'%.04f')];
                text(.35*mvf,.7*mvf,txt,'FontSize', 18)
                print('-f3','FigAllRegression.png','-dpng')
        end
end

switch element
    case 'Co'
        Conet = net;
        save Conet
    case 'Cr'
        Crnet = net;
        save Crnet
    case 'Fe'
        Fenet = net;
        save Fenet
    case 'Ni'
        Ninet = net;
        save Ninet
    case 'All'
        Allnet = net;
        save Allnet        
end

%%

% analyzeNetwork(net)


YTestPredicted = predict(net,XTest);

predictionError = YTest - YTestPredicted;

thr = 0.10;
numCorrect = sum(abs(predictionError) < thr);
numValidationImages = numel(YTest);

accuracy = numCorrect/numValidationImages

% Use the root-mean-square error (RMSE) to measure the differences between the predicted and actual angles of rotation.

squares = predictionError.^2;
rmse = sqrt(mean(squares))

% Compute R2 coefficient
Rsq = 1 - sum((YTest - YTestPredicted).^2)/sum((YTest - mean(YTest)).^2);

%

figure(30)
dx = 1.5*rmse;
plot([0 10], [0 10],'r--')
hold on
plot([0 10], [-dx 10-dx],'r--')
plot([0 10], [ dx 10+dx],'r--')
scatter((YTestPredicted'),(YTest'),50,CTest,'filled')
hold on
% scatter((YTestPredicted(:,2)'),(YTest(:,2)'),50,CTest,'filled')
xlabel("Predicted Value")
ylabel("True Value")
title("Test Set")

grid on
axis([0 1 0 1])

set(gca, 'FontSize', 18, 'LineWidth', 2)

text(.1,.75,'Formation Energy','FontSize', 18)
txt = ['RMSE ',num2str(rmse,'%0.2f')];
txt2 = ['R^2 ',num2str(Rsq,'%0.2f')];

text(.1,.7,txt,'FontSize', 18)
text(.1,.65,txt2,'FontSize', 18)

print('-f30','FigFormationTestSet.png','-dpng')