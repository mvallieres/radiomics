% *************************************************************************************
% * DESCRIPTION:                                                                      *
% * MRI study to predict different clinical outcomes in low-grade gliomas via         *   
% * texture analysis. This script computes the results presented in ref. [1].         *
% * --------------------------------------------------------------------------------- *
% * REFERENCE:                                                                        *
% * [1] Zhou, H., Vallieres, M., Bai, H.X. et al. (2017). MRI features predict        * 
% *     survival and molecular markers in diffuse lower-grade gliomas.                *
% *     Neuro-Oncology, XX(XX), 1-10. doi:10.1093/neuonc/now256                       * 
% * --------------------------------------------------------------------------------- *
% * AUTHOR(S): Martin Vallieres <mart.vallieres@gmail.com>                            *
% * --------------------------------------------------------------------------------- *
% * HISTORY:                                                                          *
% * - Creation: January 2017                                                          *
% * --------------------------------------------------------------------------------- *
% * STATEMENT:                                                                        *                                                             
% * masterScript_LGG.m: A program to predict different clinical outcomes in           * 
% * low-grade gliomas via texture analysis of MR images.                              *
% * --> Copyright (C) 2017  Martin Vallieres                                          *
% *                                                                                   *   
% *   This program is free software: you can redistribute it and/or modify            *
% *   it under the terms of the GNU General Public License as published by            *
% *   the Free Software Foundation, either version 3 of the License, or               *
% *   (at your option) any later version.                                             *
% *                                                                                   *
% *   This program is distributed in the hope that it will be useful,                 *
% *   but WITHOUT ANY WARRANTY; without even the implied warranty of                  *
% *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                   *
% *   GNU General Public License for more details.                                    *
% *                                                                                   *
% *   You should have received a copy of the GNU General Public License               *
% *   along with this program. If not, see <http://www.gnu.org/licenses/> .           *
% *************************************************************************************


% **************************** INITIALIZATIONS ****************************
clc,clear,fprintf('\n')
timeStartAll = tic;
fprintf('\n')
help masterScript_LGG
warning off
pathWORK = pwd;
pathTCGA = fullfile(pathWORK,'TCGA_DATA');
mkdir('STUDY_DATA'), pathSTUDY = fullfile(pathWORK,'STUDY_DATA');
mkdir('LOGISTIC_REGRESSION'), pathLR = fullfile(pathWORK,'LOGISTIC_REGRESSION');
mkdir('RANDOM_FORESTS'), pathRF = fullfile(pathWORK,'RANDOM_FORESTS');

% LOADING VARIABLES
cd(pathTCGA)
outcomesTCGA = load('outcomes'); outcomesTCGA = struct2cell(outcomesTCGA); outcomesTCGA = outcomesTCGA{1};
timeToEventTCGA = load('timeToEvent'); timeToEventTCGA = struct2cell(timeToEventTCGA); timeToEventTCGA = timeToEventTCGA{1};
clinicalTCGA = load('clinical'); clinicalTCGA = struct2cell(clinicalTCGA); clinicalTCGA = clinicalTCGA{1};
vasariTCGA = load('vasari'); vasariTCGA = struct2cell(vasariTCGA); vasariTCGA = vasariTCGA{1}; 
patientID = load('patientID'); patientID = struct2cell(patientID); patientID = patientID{1};
load('T1Wpath'), load('T1CEpath'), load('T2Wpath'), load('T2Fpath')
nameOutcomes = fieldnames(outcomesTCGA);
nOutcomes = length(nameOutcomes);
nPatient = numel(outcomesTCGA.(nameOutcomes{1}));

% TEXTURE EXTRACTION PARAMETERS AND DEGREES OF FREEDOM
scale_mat = [0.5,1,2,3,4]; % In mm
algo_cell = {'Equal','Uniform'};
Ng_mat = [8,16,32,64];
paramSEP = {scale_mat,algo_cell,Ng_mat}; 
baselineSEP = [2,2,3];
freedomSEP = [1,1,1];
nonTextName = {'Size','Eccentricity','Volume','Solidity'}; nNonText = numel(nonTextName);
textType = {'Global','GLSZM','GLRLM','GLCM','NGTDM'}; nTextType = numel(textType);
textName = {{'Variance','Skewness','Kurtosis'}, ...
            {'SZE','LZE','GLN','ZSN','ZP','LGZE','HGZE','SZLGE','SZHGE','LZLGE','LZHGE','GLV','ZSV'}, ...
            {'SRE','LRE','GLN','RLN','RP','LGRE','HGRE','SRLGE','SRHGE','LRLGE','LRHGE','GLV','RLV'}, ...
            {'Energy','Contrast','Entropy','Homogeneity','Correlation','SumAverage','Variance','Dissimilarity'}, ...
            {'Coarseness','Contrast','Busyness','Complexity','Strength'}};
            
            % ************************************************************* % 
            % NOTE: ADD AUTOCORRELATION AFTER THE ROI MASKS ARE UPLOADED ON %
            % THE TCIA WEBSITE.                                             % 
            % ************************************************************* %

nText = 0;
for t = 1:nTextType
    nText = nText + numel(textName{t});
end

% MULTIVARIABLE ANALYSIS PARAMETERS
nBoot = 100;
alpha = 0.5; delta = 0.5;
setSize = 25;
scans = {'T1W','T1CE','T2W','T2F'}; nScans = numel(scans);
fSetNames = {'T1W_T2W','T1W_T2F','T1CE_T2W','T1CE_T2F'}; nFset = numel(fSetNames);
maxOrder = 10;
imbalance = 'IALR';
cd(pathWORK)
nTrees = 500; % Number of random-forest trees.

% USE THESE TWO LINES (100-101) TO MANUALLY FIND THE MINE.jar FILE ON YOUR LINUX
% SYSTEM. IT COULD TAKE A LONG TIME, SO THE OTHER OPTION IS TO MANUALLY
% ENTER THE FULL PATH SIMILARLY TO WHAT YOU SEE ON LINE 102 OF THIS FILE.

tic, fprintf('\n --> INITIALIZATION: FINDING PATH TO "MINE.jar" EXECUTABLE ON THE SYSTEM ... ')
[pathMINE] = findMINE('Linux'); fprintf('DONE!\n'), toc
%pathMINE = '/home/martin/GitHub/radiomics/MultivariableModeling/MINE';

% PARALLEL COMPUTATION OPTIONS
seed =  54288; rng(seed); % For reproducibility of results. A seed chosen with heart, to remind us that it was very close and that we must continue to believe in our dream --> #1995referendum. Replace this line by "rng('shuffle')" to make it completely random.
seeds = ceil(1000000000*rand(4,1)); % A seed for feature set reduction, feature selection, prediction performance estimation and computation of final regression coefficients
nBatch = 8; % Number of parallel batch to use.
matlabPATH = 'matlab'; % Full path to the matlab executable on the system. Here, a symbolic link to the full MATLAB path has previously been created on Martin Vallieres' computer.

% DISPLAYING OPTIONS
display = false; % If this is set to true, all figures will be displayed. Setting it to false is useful to run this script in the background from terminal commands; in this case, manually open figures using "openfig('nameFigure.fig','new','visible');" after computation is over.
if ~display
    cd(pathWORK), mkdir('FIGURES'), cd('FIGURES'), pathFig = pwd; cd(pathWORK); % "FIGURES": folder where all figures will be saved
else
    pathFig = '';
end
% -------------------------------------------------------------------------




% *********************************** COMPUTATION OF RADIOMICS FEATURES FOR LGG DATA *************************************
tStart = tic;
fprintf('\n\n************************* COMPUTATION OF RADIOMICS FEATURES FOR LGG DATA *************************')

% ************ TEMPORARILY COMMENTED WHILE ROI MASKS ARE BEING UPLOADED ON THE TCIA WEBSITE ************
% % 1. READ DATA DOWNLOADED FROM THE TCIA WEBSITE (xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx)
% tic, fprintf('\n--> READING AND PROCESSING DICOM DATA FROM TCIA WEBSITE ... ')
% cd(pathTCGA), mkdir('IMAGES'), pathIMAGES = fullfile(pathTCGA,'IMAGES'); % Where we will store images and masks in MATLAB format -->  /WORKSPACE/TCGA_DATA/IMAGES
% pathDICOM = fullfile(pathTCGA,'DOI'); pathROI = fullfile(pathTCGA,'ROI'); % Where we read data from TCGA --> /WORKSPACE/TCGA_DATA/DOI for DICOM images, /WORKSPACE/TCGA_DATA/ROI for ROI
% readAllTCGA_LGG(pathDICOM,pathROI,pathIMAGES,patientID,scans,T1Wpath,T1CEpath,T2Wpath,T2Fpath)
% fprintf('DONE!\n'), toc

% ************ TEMPORARILY COMMENTED WHILE ROI MASKS ARE BEING UPLOADED ON THE TCIA WEBSITE ************
% % 2. COMPUTING NON-TEXTURE FEATURES
% tic
% cd(pathTCGA), mkdir('NON_TEXTURES'), pathNonText = fullfile(pathTCGA,'NON_TEXTURES');
% fprintf('\n--> COMPUTING NON-TEXTURE FEATURES ON %u CORES ... ',nBatch) 
% calcAllNonTextureFeatures_batchLGG(pathIMAGES,pathNonText,outcomesTCGA,nBatch,matlabPATH)
% fprintf('DONE!\n'), toc

% ************ TEMPORARILY COMMENTED WHILE ROI MASKS ARE BEING UPLOADED ON THE TCIA WEBSITE ************
% % 3. COMPUTING TEXTURE FEATURES
% tic, fprintf('\n')
% cd(pathTCGA), mkdir('TEXTURES'), pathText = fullfile(pathTCGA,'TEXTURES');
% fprintf('\n--> COMPUTING TEXTURE FEATURES ON %u CORES ... ',nBatch)
% calcAllTextures_batchLGG(pathIMAGES,pathText,outcomesTCGA,scale_mat,algo_cell,Ng_mat,nBatch,matlabPATH)
% fprintf('DONE!\n'), toc

% 4. ORGANIZING DATA
tic, fprintf('\n')
fprintf('\n--> ORGANIZING DATA FOR LGG STUDY ... ')
organizeData_LGG(pathTCGA,pathSTUDY,patientID,scans,outcomesTCGA,timeToEventTCGA,textType,textName,paramSEP);
fprintf('DONE!\n'), toc

time = toc(tStart);
fprintf('\n\nTOTAL TIME FOR COMPUTATION OF RADIOMICS FEATURES FOR LGG DATA: %f seconds\n',time)
fprintf('-------------------------------------------------------------------------------------')
% -----------------------------------------------------------------------------------------------------------



% *************** RESULTS -- PART A: BUILDING RADOMICS PREDICTION MODELS ***************
tStart = tic;
fprintf('\n\n\n********************* RESULTS -- PART A: BUILDING RADOMICS PREDICTION MODELS *********************')

% A.1 ORGANIZING TRAINING SETS
tic, fprintf('\n--> ORGANIZING TRAINING SETS  ... ')
cd(pathWORK), organizeRadiomicsExp_LGG(pathSTUDY,pathLR,nonTextName,fSetNames,paramSEP,baselineSEP,freedomSEP,textType,textName)
fprintf('DONE!\n'), toc

% A.2 PERFORM FEATURE SET REDUCTION FROM TRAINING SETS 
tic, fprintf('\n--> PERFORMING FEATURE SET REDUCTION FROM THE TRAINING SETS WITH %u CORES ... ',nBatch)
calcAllFeatureSets_batchLGG(pathMINE,pathLR,setSize,alpha,delta,nBoot,nBatch,matlabPATH,seeds(1))
fprintf('DONE!\n'), toc

% A.3 PERFORM FEATURE SET SELECTION FROM TRAINING SETS 
tic, fprintf('\n--> PERFORMING FEATURE SELECTION FROM THE TRAINING SETS WITH %u CORES ... ',nBatch)
computeAllModelChoice_batchLGG(pathLR,maxOrder,nBoot,imbalance,nBatch,matlabPATH,seeds(2))
fprintf('DONE!\n'), toc

% A.4 PERFORM PREDICTION PERFORMANCE ESTIMATION FROM TRAINING SETS 
tic, fprintf('\n--> PERFORMING PREDICTION PERFORMANCE ESTIMATION FROM THE TRAINING SETS WITH %u CORES ... ',nBatch)
computeAllPrediction_batchLGG(pathLR,maxOrder,nBoot,imbalance,nBatch,matlabPATH,seeds(3))
fprintf('DONE!\n'), toc
 
% A.5 CHOICE OF BEST PARSIMONIOUS MODELS FROM TRAINING SETS (requires user inputs)
tChoice = tic;
%chooseBestModels_LGG(pathLR,fSetNames,nameOutcomes,'AUC632',maxOrder) % USE THIS LINE TO CHOOSE YOURSELF THE BEST MODEL ORDERS
autoChoiceBestModels_LGG(pathLR,fSetNames,nameOutcomes,'AUC632',maxOrder,pathFig) % USE THIS LINE TO REPRODUCE THE RESULTS OF THIS STUDY.
timeChoice = toc(tChoice);

% A.6 COMPUTATION OF FINAL COEFFICIENTS OF RADIOMICS MODELS FROM TRAINING SETS
tic, fprintf('\n--> COMPUTING THE FINAL REGRESSION COEFFICIENTS FROM THE TRAINING SETS WITH %u CORES ... ',nBatch)
computeModelCoefficients_batchLGG(pathLR,imbalance,nBatch,matlabPATH,seeds(4))
fprintf('DONE!\n'), toc

% A.7 COMPUTATION OF FINAL COEFFICIENTS OF RADIOMICS MODELS FROM TRAINING SETS
tic, fprintf('\n--> DISPLAYING SIGMOIDAL PLOTS ... ')
displayFinalModels_LGG(pathLR,nameOutcomes,fSetNames,pathFig)
fprintf('DONE!\n'), toc

time = toc(tStart);
fprintf('\n\nTOTAL TIME FOR RESULTS -- PART A (without timeChoice): %f seconds\n',time - timeChoice)
fprintf('-------------------------------------------------------------------------------------')
% -----------------------------------------------------------------------------------------------------------



% *************** RESULTS -- PART B. COMPUTE FEATURE SELECTION FOR VASARI FEATURES **************
tStart = tic;
fprintf('\n\n\n************** RESULTS -- PART B. COMPUTE FEATURE SELECTION FOR VASARI FEATURES ***************')

% B.1 PREPARATON OF FEATURE SETS FOR VASARI FEATURES
tic, fprintf('\n--> PREPARATON OF FEATURE SETS FOR VASARI FEATURES ... ')
prepareFsetVASARI(pathWORK,pathSTUDY)
fprintf('DONE!\n'), toc

% B.2 PERFORM FEATURE SET SELECTION FOR VASARI FEATURES
tic, fprintf('\n--> PERFORM FEATURE SET SELECTION FOR VASARI FEATURES WITH %u CORES ... ',nBatch)
computeAllModelChoice_batchVASARI_LGG(fullfile(pathWORK,'VASARI'),maxOrder,nBoot,imbalance,nBatch,matlabPATH,seeds(2))
fprintf('DONE!\n'), toc

% B.3 PERFORM PREDICTION PERFORMANCE ESTIMATION FOR VASARI FEATURES
tic, fprintf('\n--> PERFORM PREDICTION PERFORMANCE ESTIMATION FOR VASARI FEATURES WITH %u CORES ... ',nBatch)
computeAllPrediction_batchVASARI_LGG(fullfile(pathWORK,'VASARI'),maxOrder,nBoot,imbalance,nBatch,matlabPATH,seeds(3))
fprintf('DONE!\n'), toc

% B.4 CHOICE OF BEST PARSIMONIOUS MODELS FOR VASARI FEATURES
tChoice = tic;
%chooseBestModels_LGG(fullfile(pathWORK,'VASARI'),{'VASARI'},nameOutcomes,'AUC632',maxOrder) % USE THIS LINE TO CHOOSE YOURSELF THE BEST MODEL ORDERS
autoChoiceBestModels_VASARI_LGG(fullfile(pathWORK,'VASARI'),{'VASARI'},nameOutcomes,'AUC632',maxOrder,pathFig) % USE THIS LINE TO REPRODUCE THE RESULTS OF THIS STUDY.
timeChoice = toc(tChoice);

% B.5 PRODUCE SINGLE PLOT FOR VASARI FEATURES
tic, fprintf('\n--> DISPLAYING VASARI PLOT ... ')
produceVASARIplot_LGG(fullfile(pathWORK,'VASARI','RESULTS'),nameOutcomes,maxOrder,pathFig)
fprintf('DONE!\n'), toc

time = toc(tStart);
fprintf('\n\nTOTAL TIME FOR RESULTS -- PART B: %f seconds\n',time - timeChoice)
fprintf('-------------------------------------------------------------------------------------')
% -------------------------------------------------------------------------



% *************** RESULTS -- PART C. COMPUTE RANDOM FORESTS EXPERIMENTS **************
tStart = tic;
fprintf('\n\n\n************** RESULTS -- PART C. COMPUTE RANDOM FORESTS EXPERIMENTS ***************')

% C.1 ORGANIZING RANDOM FORESTS EXPERIMENTS
tic, fprintf('\n--> ORGANIZING RANDOM FORESTS EXPERIMENTS ... ')
organizeRF_LGG(pathTCGA,pathSTUDY,fullfile(pathWORK,'LOGISTIC_REGRESSION','FINAL_MODELS'),fullfile(pathWORK,'VASARI','FINAL_MODELS'),pathRF);
fprintf('DONE!\n'), toc

% C.2 PERFORM RANDOM FORESTS EXPERIMENTS: COMBINED (TEXTURE + VASARI + CLINICAL)
tic, fprintf('\n--> PERFORM RANDOM FORESTS EXPERIMENTS: TEXTURE + VASARI + CLINICAL ... ')
performBatchRF_LGG(pathRF,nTrees,seeds(3),matlabPATH)
fprintf('DONE!\n'), toc

% C.3 PERFORM RANDOM FORESTS EXPERIMENTS: IMAGING (TEXT + VASARI) ONLY
tic, fprintf('\n--> PERFORM RANDOM FORESTS EXPERIMENTS: TEXTURE + VASARI ... ')
performBatchRFimaging_LGG(pathRF,nTrees,seeds(3),matlabPATH)
fprintf('DONE!\n'), toc

% C.4 PERFORM RANDOM FORESTS EXPERIMENTS: CLINICAL ONLY
tic, fprintf('\n--> PERFORM RANDOM FORESTS EXPERIMENTS: CLINICAL ONLY ... ')
performBatchRFclinical_LGG(pathRF,nTrees,seeds(3),matlabPATH)
fprintf('DONE!\n'), toc

% C.5 PERFORM RANDOM FORESTS EXPERIMENTS: TEXTURES ONLY
tic, fprintf('\n--> PERFORM RANDOM FORESTS EXPERIMENTS: TEXTURES ONLY ... ')
performBatchRFtextures_LGG(pathRF,nTrees,seeds(3),matlabPATH)
fprintf('DONE!\n'), toc

% C.6 PERFORM RANDOM FORESTS EXPERIMENTS: VASARI ONLY
tic, fprintf('\n--> PERFORM RANDOM FORESTS EXPERIMENTS: VASARI ONLY ... ')
performBatchRFvasari_LGG(pathRF,nTrees,seeds(3),matlabPATH)
fprintf('DONE!\n'), toc

% C.7 PERFORM RANDOM FORESTS EXPERIMENTS: TEXTURE + CLINICAL
tic, fprintf('\n--> PERFORM RANDOM FORESTS EXPERIMENTS: TEXTURES + CLINICAL ... ')
performBatchRFtextClin_LGG(pathRF,nTrees,seeds(3),matlabPATH)
fprintf('DONE!\n'), toc

time = toc(tStart);
fprintf('\n\nTOTAL TIME FOR RESULTS -- PART C: %f seconds\n',time)
fprintf('-------------------------------------------------------------------------------------')
% -------------------------------------------------------------------------



timeAll = toc(timeStartAll);
fprintf('\n\n\n---------> TOTAL TIME FOR THIS STUDY: %f hours',timeAll/3600)
fprintf('\n\n\n**************************************************************************\n')
fprintf('------------------------------- THE END ----------------------------------\n')
fprintf('**************************************************************************\n\n')
% **************************************************************************
% ------------------------------- THE END ----------------------------------
% **************************************************************************
