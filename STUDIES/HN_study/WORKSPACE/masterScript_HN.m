% *************************************************************************************
% * DESCRIPTION:                                                                      *
% * This script computes all the experiments presented in ref. [1].                   *
% * --------------------------------------------------------------------------------- *
% * AUTHOR(S): Martin Vallieres <mart.vallieres@gmail.com>                            *
% * --------------------------------------------------------------------------------- *
% * HISTORY:                                                                          *
% * - Creation: August 2016                                                           *
% * --------------------------------------------------------------------------------- *
% * STATEMENT:                                                                        *
% * This file is part of <https://github.com/mvallieres/radiomics/>,                  *
% * a package providing MATLAB programming tools for radiomics analysis.              *
% * --> Copyright (C) 2015, 2016  Martin Vallieres                                    *
% *                                                                                   *
% *   This package is free software: you can redistribute it and/or modify            * 
% *   it under the terms of the GNU General Public License as published by            *
% *   the Free Software Foundation, either version 3 of the License, or               *
% *   (at your option) any later version.                                             *
% *                                                                                   *
% *   This package is distributed in the hope that it will be useful,                 *
% *   but WITHOUT ANY WARRANTY; without even the implied warranty of                  *
% *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                   *
% *   GNU General Public License for more details.                                    *
% *                                                                                   *
% *   You should have received a copy of the GNU General Public License               *
% *   along with this package.  If not, see <http://www.gnu.org/licenses/>.           *
% *************************************************************************************

%%%%%%%%%%%%%%%%%%%%% LIST OF IMPORTANT ABBREVIATIONS %%%%%%%%%%%%%%%%%%%%%
% - LR: Logistic regression                                               %
% - CR: Cox regression                                                    %
% - RF: Random forests                                                    %
% - sign: Radiomics signature of (Aerts et al., Nat. Commun., 2014)       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% **************************** INITIALIZATIONS ****************************
clc,clear,fprintf('\n')
timeStartAll = tic;
help masterScript_HN
warning off
cohorts = {'HGJ','HMR','CHUS','CHUM'}; nCohort = numel(cohorts);
featType = {'GTVp','GTVtot'}; nFeatType = numel(featType);

% LOADING VARIABLES
pathWORK = pwd;
clinical = load('clinical'); clinical = struct2cell(clinical); clinical = clinical{1}; % Clinical parameters: Age, T, N, TNM, HPV
load('outcomes'), load('roiNames'), load('timeToEvent'), load('subTypes') % Variables 'outcomes' and 'roiNames' now in the workspace
count = 0;
for i = 1:nCohort
    nPatient.(cohorts{i}) = size(roiNames.(cohorts{i}),1);
    count = count + nPatient.(cohorts{i});
end
nPatient.TOTAL = count;
nameOutcomes = fieldnames(outcomes.(cohorts{1})); 
nOutcomes = numel(nameOutcomes);

% TEXTURE EXTRACTION PARAMETERS AND DEGREES OF FREEDOM
scale_mat = [1,2,3,4,5];
algo_cell = {'Equal','Uniform'};
Ng_mat = [8,16,32,64];
paramSEP = {scale_mat,algo_cell,Ng_mat};
nameSEP = {'Scale','Quant.algo','Ng'};
baselineSEP = [3,2,3];
freedomSEP = [1,1,1];
nonTextName = {'SUVmax','SUVpeak','SUVmean','aucCSH','TLG','PercentInactive','gETU','Volume','Size','Solidity','Eccentricity','Compactness'}; nNonText = numel(nonTextName);
textType = {'Global','GLCM','GLRLM','GLSZM','NGTDM'}; nTextType = numel(textType);
textName = {{'Variance','Skewness','Kurtosis'}, ...
            {'Energy','Contrast','Entropy','Homogeneity','Correlation','SumAverage','Variance','Dissimilarity','AutoCorrelation'}, ...
            {'SRE','LRE','GLN','RLN','RP','LGRE','HGRE','SRLGE','SRHGE','LRLGE','LRHGE','GLV','RLV'}, ...
            {'SZE','LZE','GLN','ZSN','ZP','LGZE','HGZE','SZLGE','SZHGE','LZLGE','LZHGE','GLV','ZSV'}, ...
            {'Coarseness','Contrast','Busyness','Complexity','Strength'}};
nText = 0;
for t = 1:nTextType
    nText = nText + numel(textName{t});
end
        
% MULTIVARIABLE ANALYSIS PARAMETERS
nBoot = 100; % Number of bootstrapping experiments used in models construction
alpha = 0.5; delta = 0.5; setSize = 25; % Feature set reduction parameters. See (Vallières et al., Phys. Med. Biol., 2015) for more details
fSetNames = {'PET','CT','PETCT'}; nFset = numel(fSetNames); % PETCT --> Separate PET and CT radiomic features combined into one set
maxOrder = 10; % Maximum number or radiomics variable combinations
imbalance = 'IALR'; % Imbalance-adjusted logistic regression
tic, fprintf('\n --> INITIALIZATION: FINDING PATH TO "MINE.jar" EXECUTABLE ON THE SYSTEM ... ')
[pathMINE] = findMINE('Linux'); fprintf('DONE!\n'), toc
testCost = 0.5:0.1:2; % Emphasis factor on positive instances during random forest training
testSplit = 1/3; % Proportion of test cases in stratified random sub-sampling splits

% PARALLEL AND RANDOMIZATION OPTIONS
seed =  54288; rng(seed); % For reproducibility of results. A seed chosen with heart, to remind us that it was very close and that we must continue to believe in our dream --> #1995referendum. Replace this line by "rng('shuffle')" to make it completely random. 
seeds = ceil(1000000000*rand(4,1)); % A bootstrapping seed for feature set reduction, feature selection, prediction performance estimation and computation of final regression coefficients
nBatch = 9; % Number of parallel batch to use (3 outcomes * 3 different feature sets)
nBatch_Read = 4; % beware: RAM usage limitations
matlabPATH = 'matlab'; % Full path to the matlab executable on the system. Here, a symbolic link to the full MATLAB path has previously been created on Martin Vallieres' computer.

% DISPLAYING OPTIONS
display = false; % If this is set to true, all figures will be displayed. Setting it to false is useful to run this script in the background from terminal commands; in this case, manually open figures using "openfig('nameFigure.fig','new','visible');" after computation is over.
if ~display
    cd(pathWORK), mkdir('FIGURES'), cd('FIGURES'), pathFig = pwd; cd(pathWORK); % "FIGURES": folder where all figures will be saved
else
    pathFig = '';
end
% -------------------------------------------------------------------------



% *********************************** COMPUTATION OF RADIOMICS FEATURES *************************************
tStart = tic;
fprintf('\n\n************************* COMPUTATION OF RADIOMICS FEATURES *************************')
          
% 1. READ DATA DOWNLOADED FROM THE TCIA WEBSITE (xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx)
tic, fprintf('\n--> READING AND PROCESSING DICOM DATA FROM TCIA WEBSITE ON %u CORES ... ',nBatch_Read)
mkdir('DATA'), pathData = fullfile(pathWORK,'DATA'); cd(fullfile(pathWORK,'DICOM'))
readAllDICOM(fullfile(pathWORK,'DICOM'),pathData,nBatch_Read)
fprintf('DONE!\n'), toc

% 2. ORGANIZING PATIENT DATA BY PET, CT and CTsim NAMES
tic, fprintf('\n--> ORGANIZING PATIENT DATA BY PET, CT and CTsim NAMES ... ')
for i = 1:nCohort
    correctPatientNames_HN(pathData,cohorts{i},1:(nPatient.(cohorts{i})))
    patientNamePT.((cohorts{i})) = cell(nPatient.(cohorts{i}),1);
    patientNameCT.((cohorts{i})) = cell(nPatient.(cohorts{i}),1);
    for j = 1:nPatient.(cohorts{i})
        patientNamePT.((cohorts{i})){j} = [cohorts{i},'_',num2str(j,'%.3i'),'_PT.PTscan.mat'];
        patientNameCT.((cohorts{i})){j} = [cohorts{i},'_',num2str(j,'%.3i'),'_CT.CTscan.mat'];
    end
end
fprintf('DONE!\n'), toc
 
% 3. COMPUTING NON-TEXTURE FEATURES
tic
cd(pathWORK), mkdir('FEATURES'), pathFeatures = fullfile(pathWORK,'FEATURES');
cd(pathFeatures), mkdir('NON_TEXTURES'), pathNonText = fullfile(pathFeatures,'NON_TEXTURES');
for i = 1:nCohort
   cohort = cohorts{i}; namePT = patientNamePT.(cohort); nameCT = patientNameCT.(cohort); outcomeCohort = outcomes.(cohort);
   fprintf('\n--> COMPUTING NON-TEXTURE FEATURES OF ''%s'' COHORT ON %u CORES ... ',cohort,nBatch)
   for j = 1:nFeatType
       nameROI = roiNames.(cohort)(:,j); type = featType{j}; 
       calcAllNonTextureFeatures_batchHN(pathData,pathNonText,namePT,nameCT,nameROI,outcomeCohort,type,nBatch,matlabPATH)
   end
   fprintf('DONE!')
end
fprintf('\n'), toc

% 4. COMPUTING PET AND CT TEXTURE FEATURES
tic, fprintf('\n')
cd(pathFeatures), mkdir('TEXTURES'), pathText = fullfile(pathFeatures,'TEXTURES');
for i = 1:nCohort
   cohort = cohorts{i}; namePT = patientNamePT.(cohort); nameCT = patientNameCT.(cohort); outcomeCohort = outcomes.(cohort);
   fprintf('\n--> COMPUTING TEXTURE FEATURES OF PET AND CT SCANS OF "%s" COHORT ON %u CORES ... ',cohort,nBatch)
   for j = 1:nFeatType
       nameROI = roiNames.(cohort)(:,j); type = featType{j};
       calcAllSeparateTextures_batchHN(pathData,pathText,namePT,nameCT,nameROI,outcomeCohort,type,scale_mat,algo_cell,Ng_mat,nBatch,matlabPATH)
   end
   fprintf('DONE!')
end
fprintf('\n'), toc

time = toc(tStart);
fprintf('\n\nTOTAL TIME FOR COMPUTATION OF RADIOMICS FEATURES: %f seconds\n',time)
fprintf('-------------------------------------------------------------------------------------')
% -----------------------------------------------------------------------------------------------------------



% *************** RESULTS -- PART A: ASSOCIATION OF RADIOMICS AND CLINICAL VARIABLES WITH OUTCOMES ***************
tStart = tic;
fprintf('\n\n\n********************* RESULTS -- PART A: ASSOCIATION OF RADIOMICS AND CLINICAL VARIABLES WITH OUTCOMES *********************')
cd(pathWORK), mkdir('UNIVARIATE_RESULTS'), cd('UNIVARIATE_RESULTS'), pathUnivariate = pwd; cd(pathWORK)

% A.1 COMPUTING ASSOCIATIONS BETWEEN CLINICAL VARIABLES AND OUTCOMES
tic, fprintf('\n--> COMPUTING ASSOCIATIONS BETWEEN CLINICAL VARIABLES AND OUTCOMES ... ')
computeAssociations_Clinical(pathUnivariate,cohorts,clinical,outcomes) 
fprintf('DONE!\n'), toc

% A.2 COMPUTING ASSOCIATIONS BETWEEN NON-TEXTURE FEATURES AND OUTCOMES
tic, fprintf('\n--> COMPUTING ASSOCIATIONS BETWEEN NON-TEXTURE FEATURES AND OUTCOMES ... ')
computeAssociations_NonTextures(pathUnivariate,pathFeatures,cohorts,outcomes,featType,nonTextName)
fprintf('DONE!\n'), toc

% A.3 COMPUTING ASSOCIATIONS BETWEEN TEXTURE FEATURES AND OUTCOMES
tic, fprintf('\n--> COMPUTING ASSOCIATIONS BETWEEN TEXTURE FEATURES AND OUTCOMES ... ')
computeAssociations_Textures(pathUnivariate,pathFeatures,cohorts,outcomes,featType,textType,textName,nPatient,scale_mat,algo_cell,Ng_mat)
fprintf('DONE!\n'), toc

time = toc(tStart);
fprintf('\n\nTOTAL TIME FOR PART RESULTS -- PART A: %f seconds\n',time)
fprintf('-------------------------------------------------------------------------------------')
% -----------------------------------------------------------------------------------------------------------
%%% FROM HERE, ONLY THE "GTVtot" ROI IS USED %%%



% *************** RESULTS -- PART B: BUILDING RADOMICS PREDICTION MODELS ***************
tStart = tic;
fprintf('\n\n\n********************* RESULTS -- PART B: BUILDING RADOMICS PREDICTION MODELS *********************')

% B.1 ORGANIZING TRAINING ("HGJ + CHUS") AND TESTING ("HMR + CHUM") SETS 
tic, fprintf('\n--> ORGANIZING TRAINING ("HGJ + CHUS") AND TESTING ("HMR + CHUM") SETS  ... ')
cd(pathWORK), mkdir('LOGISTIC_REGRESSION'), cd('LOGISTIC_REGRESSION'), pathLR = pwd; cd(pathWORK), mkdir('COX_REGRESSION'), cd('COX_REGRESSION'), pathCR = pwd;
cd(pathWORK), organizeExperiments_Radiomics(pathFeatures,pathLR,pathCR,outcomes,timeToEvent,nPatient,nonTextName,fSetNames,paramSEP,baselineSEP,freedomSEP,textType,textName)
fprintf('DONE!\n'), toc

% B.2 PERFORM FEATURE SET REDUCTION FROM TRAINING SETS 
tic, fprintf('\n--> PERFORMING FEATURE SET REDUCTION FROM THE TRAINING SETS WITH %u CORES ... ',nBatch)
calcAllFeatureSets_batchHN(pathMINE,pathLR,setSize,alpha,delta,nBoot,nBatch,matlabPATH,seeds(1))
fprintf('DONE!\n'), toc

% B.3 PERFORM FEATURE SET SELECTION FROM TRAINING SETS 
tic, fprintf('\n--> PERFORMING FEATURE SELECTION FROM THE TRAINING SETS WITH %u CORES ... ',nBatch)
computeAllModelChoice_batchHN(pathLR,maxOrder,nBoot,imbalance,nBatch,matlabPATH,seeds(2))
fprintf('DONE!\n'), toc

% B.4 PERFORM PREDICTION PERFORMANCE ESTIMATION FROM TRAINING SETS 
tic, fprintf('\n--> PERFORMING PREDICTION PERFORMANCE ESTIMATION FROM THE TRAINING SETS WITH %u CORES ... ',nBatch)
computeAllPrediction_batchHN(pathLR,maxOrder,nBoot,imbalance,nBatch,matlabPATH,seeds(3))
fprintf('DONE!\n'), toc
 
% B.5 CHOICE OF BEST PARSIMONIOUS MODELS FROM TRAINING SETS (requires user inputs)
tChoice = tic;
%chooseBestModels_HN(pathLR,pathCR,fSetNames,nameOutcomes,'AUC632',maxOrder) % USE THIS LINE TO CHOOSE YOURSELF THE BEST MODEL ORDERS
autoChoiceBestModels_HN(pathLR,pathCR,fSetNames,nameOutcomes,'AUC632',maxOrder,pathFig) % USE THIS LINE TO REPRODUCE THE RESULTS OF THIS STUDY.
timeChoice = toc(tChoice);

% B.6 COMPUTATION OF FINAL COEFFICIENTS OF RADIOMICS MODELS FROM TRAINING SETS
tic, fprintf('\n--> COMPUTING THE FINAL REGRESSION COEFFICIENTS FROM THE TRAINING SETS WITH %u CORES ... ',nBatch)
computeModelCoefficients_batchHN(pathLR,imbalance,nBatch,matlabPATH,seeds(4))
computeModelCoefficientsTime_batchHN(pathCR,nBatch,matlabPATH,seeds(4))
fprintf('DONE!\n'), toc
 
time = toc(tStart);
fprintf('\n\nTOTAL TIME FOR RESULTS -- PART B (without timeChoice): %f seconds\n',time - timeChoice)
fprintf('-------------------------------------------------------------------------------------')
% -----------------------------------------------------------------------------------------------------------
 
 
 
% *************** RESULTS -- PART C: PROGNOSTIC TESTING OF RADIOMICS MODELS ***************
tStart = tic;
fprintf('\n\n\n********************* RESULTS -- PART C: PROGNOSTIC TESTING OF RADIOMICS MODELS *********************')

% C.1 TESTING THE RADIOMICS MODELS
tic, fprintf('\n--> TESTING THE RADIOMICS MODELS ...')
testFinalModels_HN(pathLR,fSetNames,{paramSEP,paramSEP,paramSEP},{nameSEP,nameSEP,nameSEP})
testFinalModelsTime_HN(pathCR,fSetNames,{paramSEP,paramSEP,paramSEP},{nameSEP,nameSEP,nameSEP})
fprintf('DONE!\n'), toc

% C.2 DISPLAYING THE PREDICTION PERFORMANCE OF RADIOMICS MODELS
tic, fprintf('\n--> DISPLAYING THE PREDICTION PERFORMANCE OF RADIOMICS MODELS ...')
for o = 1:nOutcomes
    nameOutcome = nameOutcomes{o};
    displayBarResults_MetricBased_HN(pathLR,fSetNames,nameOutcome,{'AUC','Sensitivity','Specificity','Accuracy'},'LR',pathFig)
end
fprintf('DONE!\n'), toc

% C.3 DISPLAYING STANDARD KAPLAN-MEIER CURVES FOR RADIOMIC MODELS
tic, fprintf('\n--> DISPLAYING STANDARD KAPLAN-MEIER CURVES FOR RADIOMIC MODELS ...')
displayKaplanMeierCR_HN(pathCR,nameOutcomes,{'PET','PET','PET'},pathFig)
displayKaplanMeierCR_HN(pathCR,nameOutcomes,{'CT','CT','CT'},pathFig)
displayKaplanMeierCR_HN(pathCR,nameOutcomes,{'PETCT','PETCT','PETCT'},pathFig)
fprintf('DONE!\n'), toc

time = toc(tStart);
fprintf('\n\nTOTAL TIME FOR RESULTS -- PART C: %f seconds\n',time)
fprintf('-------------------------------------------------------------------------------------')
% -----------------------------------------------------------------------------------------------------------



% *************** RESULTS -- PART D: COMBINATION OF RADIOMICS MODELS WITH CLINICAL VARIABLES ***************
tStart = tic;
fprintf('\n\n\n************** RESULTS -- PART D: COMBINATION OF RADIOMICS MODELS WITH CLINICAL VARIABLES ***************')
cd(pathWORK), mkdir('RANDOM_FORESTS'), cd('RANDOM_FORESTS'), pathRF = pwd; cd(pathWORK)

% D.1 ORGANIZING RANDOM FOREST DATA
tic, fprintf('\n--> ORGANIZING RANDOM FOREST DATA ... ')
organizeRFexperiments_HN(pathWORK,pathLR,pathRF,nameOutcomes,fSetNames)
fprintf('DONE!\n'), toc

% D.2 ESTIMATING AND TESTING THE BEST RANDOM FORESTS "CLINICAL" VARIABLES
tic, fprintf('\n--> ESTIMATING AND AND TESTING THE BEST "CLINICAL" VARIABLES FOR RANDOM FORESTS) ...')
nSplit = 10; computeClinicalResults_RF(pathRF,nBoot,nSplit,testSplit,testCost,[seeds(3),seeds(4)],matlabPATH)
fprintf('DONE!\n'), toc

% D.3 ESTIMATING THE BEST "RADIOMICS + CLINICAL" RANDOM FOREST PARAMETERS FROM TRAINING SETS
tic, fprintf('\n--> ESTIMATING THE BEST "RADIOMICS + CLINICAL" RANDOM FOREST PARAMETERS FROM TRAINING SETS ... ')
nSplit = 10; estimateCombinedRF_batchHN(pathRF,nSplit,nBoot,testSplit,testCost,seeds(3),matlabPATH)
fprintf('DONE!\n'), toc

% D.4 COMPUTING FINAL RANDOM FORESTS CLASSIFIERS FROM TRAINING SETS
tic, fprintf('\n--> COMPUTING FINAL RANDOM FORESTS CLASSIFIERS FROM TRAINING SETS ... ')
computeAllFinalRF_HN(pathRF,nBoot,seeds(4))
fprintf('DONE!\n'), toc

% D.5 TESTING THE FINAL RANDOM FORESTS CLASSIFIERS
tic, fprintf('\n--> TESTING THE FINAL RANDOM FORESTS CLASSIFIERS ... ')
predictAllRF_HN(pathRF)
fprintf('DONE!\n'), toc

% D.6 CALCULATING THE SIGNIFICANCE OF PREDICTION IMPROVEMENTS BY ADDING CLINICAL VARIABLES
tic, fprintf('\n--> CALCULATING THE SIGNIFICANCE OF PREDICTION IMPROVEMENTS BY ADDING CLINICAL VARIABLES ... ')
calcAll_AUCcomp(pathLR,pathRF,nameOutcomes,fSetNames)
fprintf('DONE!\n'), toc

% D.7 DISPLAYING THE PREDICTION PERFORMANCE OF FINAL RANDOM FORESTS
tic, fprintf('\n--> DISPLAYING THE PREDICTION PERFORMANCE OF FINAL RANDOM FORESTS ...')
fSetNamesClin = {'PETclinic','CTclinic','PETCTclinic'};
for o = 1:nOutcomes
    nameOutcome = nameOutcomes{o};
    displayBarResults_MetricBased_HN(pathRF,fSetNamesClin,nameOutcome,{'AUC','Sensitivity','Specificity','Accuracy'},'RF',pathFig)
end
fprintf('DONE!\n'), toc

time = toc(tStart);
fprintf('\n\nTOTAL TIME FOR RESULTS -- PART D: %f seconds\n',time)
fprintf('-------------------------------------------------------------------------------------')
% -------------------------------------------------------------------------



% **************** RESULTS -- PART E: COMPUTING RESULTS FOR THE "VOLUME" VARIABLE ****************
tStart = tic;
fprintf('\n\n\n************** RESULTS -- PART E: COMPUTING RESULTS FOR THE "VOLUME" VARIABLE ***************')

% E.1 COMPUTING AND TESTING REGRESSION MODELS FOR "VOLUME" VARIABLE(S)
tic, fprintf('\n--> COMPUTING AND TESTING REGRESSION CLASSIFERS FOR "VOLUME" VARIABLE(S) ...')
computeVolumeResults_Regression(pathLR,pathCR,imbalance,seeds(4),pathFig)
fprintf('DONE!\n'), toc

% E.2 COMPUTING AND TESTING RANDOM FORESTS FOR VOLUME + CLINICAL VARIABLES
tic, fprintf('\n--> COMPUTING AND TESTING RANDOM FORESTS FOR VOLUME + CLINICAL VARIABLES ...')
cd(pathLR), load('training'), volumeTrain = training.(nameOutcomes{1}).nonText.Volume.Data; load('testing'), volumeTest = testing.(nameOutcomes{1}).nonText.Volume.Data;
cd(pathRF), save('volumeTrain','volumeTrain'); save('volumeTest','volumeTest'); cd(pathWORK)
nSplit = 10; computeVolumeResults_RF(pathRF,nBoot,nSplit,testSplit,testCost,[seeds(3),seeds(4)],matlabPATH)
fprintf('DONE!\n'), toc

time = toc(tStart);
fprintf('\n\nTOTAL TIME FOR RESULTS -- PART E: %f seconds\n',time)
fprintf('-------------------------------------------------------------------------------------')
% -------------------------------------------------------------------------



% **************** RESULTS -- PART F: COMPUTING RESULTS FOR (Aerts,2014) RADIOMIC SIGNATURE ****************
tStart = tic;
fprintf('\n\n\n************** RESULTS -- PART F: COMPUTING RESULTS FOR (Aerts,2014) RADIOMIC SIGNATURE ***************')

% F.1 COMPUTING FEATURES OF (Aerts et al. 2014) SIGNATURE
tic
cd(pathFeatures), mkdir('Aerts2014'), pathAerts = fullfile(pathFeatures,'Aerts2014');
for i = 1:nCohort
    cohort = cohorts{i}; namePT = patientNamePT.(cohort); nameCT = patientNameCT.(cohort); outcomeCohort = outcomes.(cohort);
    fprintf('\n--> COMPUTING (Aerts,2014) FEATURES ON PET AND CT SCANS OF "%s" COHORT WITH %u CORES ... ',cohort,nBatch)
    for j = 1:nFeatType
        nameROI = roiNames.(cohort)(:,j); type = featType{j};
        calcAllAertsFeatures_batchHN(pathData,pathAerts,namePT,nameCT,nameROI,outcomeCohort,type,nBatch_Read,matlabPATH)
    end
    fprintf('DONE!')
end
fprintf('\n'), toc

% F.2 ORGANIZING TRAINING ("HGJ + CHUS") AND TESTING ("HMR + CHUM") SETS 
tic, fprintf('\n--> ORGANIZING TRAINING ("HGJ + CHUS") AND TESTING ("HMR + CHUM") SETS FOR (Aerts et al. 2014) FEATURES  ... ')
cd(pathWORK), organizeExperiments_RadiomicsSign(pathFeatures,pathLR,pathCR,pathRF,outcomes)
fprintf('DONE!\n'), toc

% F.3 COMPUTATION OF FINAL COEFFICIENTS OF RADIOMICS SIGNATURE (Aerts et al. 2014) FROM TRAINING SETS
tic, fprintf('\n--> COMPUTATION OF FINAL COEFFICIENTS OF RADIOMICS SIGNATURE (Aerts et al. 2014) FROM TRAINING SETS  ... ')
computeModelCoefficientsAerts_batchHN(pathLR,imbalance,nBatch,matlabPATH,seeds(4)) % Computing new logistic regression coefficients for radiomic signature
computeModelCoefficientsTimeAerts_batchHN(pathCR,nBatch,matlabPATH,seeds(4)) % Computing new cox regression coefficients for radiomic signature
fprintf('DONE!\n'), toc

% F.4 TESTING THE RADIOMICS SIGNATURE OF (Aerts et al., 2014)
tic, fprintf('\n--> TESTING THE RADIOMICS SIGNATURE OF (Aerts et al., 2014)  ... ')
testFinalModelsAerts_HN(pathLR,{'PET','CT','CTorig'}) % Logistic regression results, with coefficients retrained in training set
testFinalModelsTimeAerts_HN(pathCR,{'PET','CT','CTorig'}) % Cox regression results, with coefficients retrained in training set
testFinalModelsTimeAertsComplete_HN(pathCR) % Cox regression results, using original coefficients from (Aerts et al., 2014)
fprintf('DONE!\n'), toc

% F.5 DISPLAYING THE PREDICTION PERFORMANCE OF RADIOMICS SIGNATURE OF (Aerts et al., 2014)
tic, fprintf('\n--> DISPLAYING THE PREDICTION PERFORMANCE OF RADIOMICS SIGNATURE OF (Aerts et al., 2014) ...')
for o = 1:nOutcomes
    nameOutcome = nameOutcomes{o};
    displayBarResults_MetricBased_HN(pathLR,{'PETsign','CTsign','CTorigsign'},nameOutcome,{'AUC','Sensitivity','Specificity','Accuracy'},'LR',pathFig)
end
fprintf('DONE!\n'), toc

% F.6 DISPLAYING STANDARD KAPLAN-MEIER CURVES OF RADIOMICS SIGNATURE OF (Aerts et al., 2014)
tic, fprintf('\n--> DISPLAYING STANDARD KAPLAN-MEIER CURVES OF RADIOMICS SIGNATURE OF (Aerts et al., 2014) ...')
displayKaplanMeierCR_HN(pathCR,nameOutcomes,{'PETsign','PETsign','PETsign'},pathFig)
displayKaplanMeierCR_HN(pathCR,nameOutcomes,{'CTsign','CTsign','CTsign'},pathFig)
displayKaplanMeierCR_HN(pathCR,nameOutcomes,{'CTorigsign','CTorigsign','CTorigsign'},pathFig)
displayKaplanMeierCR_HN(pathCR,{'Death'},{'CTorigsignComplete'},pathFig)
fprintf('DONE!\n'), toc

% F.7 ESTIMATING THE BEST "RADIOMICS + CLINICAL" RANDOM FOREST PARAMETERS FROM TRAINING SETS FOR RADIOMICS SIGNATURE OF (Aerts et al., 2014)
tic, fprintf('\n--> ESTIMATING THE BEST "RADIOMICS + CLINICAL" RANDOM FOREST PARAMETERS FROM TRAINING SETS FOR RADIOMICS SIGNATURE OF (Aerts et al., 2014)  ... ')
nSplit = 10; estimateCombinedRF_Aerts_batchHN(pathRF,nSplit,nBoot,testSplit,testCost,seeds(3),matlabPATH)
fprintf('DONE!\n'), toc

% F.8 COMPUTING FINAL RANDOM FORESTS CLASSIFIERS FROM TRAINING SETS FOR RADIOMICS SIGNATURE OF (Aerts et al., 2014)
tic, fprintf('\n--> COMPUTING FINAL RANDOM FORESTS CLASSIFIERS FROM TRAINING SETS FOR RADIOMICS SIGNATURE OF (Aerts et al., 2014) ... ')
computeAllFinalRF_Aerts(pathRF,nBoot,seeds(4))
fprintf('DONE!\n'), toc

% F.9 TESTING THE FINAL RANDOM FORESTS CLASSIFIERS FOR RADIOMICS SIGNATURE OF (Aerts et al., 2014)
tic, fprintf('\n--> TESTING THE FINAL RANDOM FORESTS CLASSIFIERS FOR RADIOMICS SIGNATURE OF (Aerts et al., 2014) ... ')
predictAllRF_Aerts(pathRF)
fprintf('DONE!\n'), toc

% F.10 DISPLAYING THE PREDICTION PERFORMANCE OF FINAL RANDOM FORESTS (Rad signature + clinical)
tic, fprintf('\n--> COMPARING THE PREDICTION PERFORMANCE OF FINAL RANDOM FORESTS (Rad signature + clinical) ...')
fSetNamesSignClinic = {'PETsignClinic','CTsignClinic','CTorigsignClinic'};
for o = 1:nOutcomes
    nameOutcome = nameOutcomes{o};
    displayBarResults_MetricBased_HN(pathRF,fSetNamesSignClinic,nameOutcome,{'AUC','Sensitivity','Specificity','Accuracy'},'RF',pathFig)
end
fprintf('DONE!\n'), toc

time = toc(tStart);
fprintf('\n\nTOTAL TIME FOR RESULTS -- PART F: %f seconds\n',time)
fprintf('-------------------------------------------------------------------------------------')
% -------------------------------------------------------------------------



% *************** RESULTS -- PART G: RISK ASSESSMENT ***************
tStart = tic;
fprintf('\n\n\n************** RESULTS -- PART G: RISK ASSESSMENT ***************')

% G.1 DISPLAYING RANDOM FORESTS OUTPUT PROBABILITY OF OCCURRENCE
tic, fprintf('\n--> DISPLAYING THE PREDICTION PERFORMANCE OF FINAL RANDOM FORESTS ...')
displayProbGraph_HN(pathRF,{'Locoregional','Distant','Death'},{'PETCTclinic','CTclinic','clinic'},[33,66],pathFig) % Only showing probability results of these three best situations: 1) PET_Locoregional, 2) CT_Distant, 3) PETCT_Death
fprintf('DONE!\n'), toc

% G.2 DISPLAYING KAPLAN-MEIER CURVES FOR 2 RISK GROUPS (STANDARD: p < 0.5, p > 0.5)
tic, fprintf('\n--> DISPLAYING KAPLAN-MEIER CURVES FOR 2 RISK GROUPS (STANDARD: p < 0.5, p > 0.5) ...')
displayKaplanMeierRF_HN(pathRF,nameOutcomes,{'clinic','clinic','clinic'},pathFig)
displayKaplanMeierRF_HN(pathRF,nameOutcomes,{'VOLclinic','VOLclinic','VOLclinic'},pathFig)
displayKaplanMeierRF_HN(pathRF,nameOutcomes,{'PETclinic','PETclinic','PETclinic'},pathFig)
displayKaplanMeierRF_HN(pathRF,nameOutcomes,{'CTclinic','CTclinic','CTclinic'},pathFig)
displayKaplanMeierRF_HN(pathRF,nameOutcomes,{'PETCTclinic','PETCTclinic','PETCTclinic'},pathFig)
displayKaplanMeierRF_HN(pathRF,nameOutcomes,{'PETsignClinic','PETsignClinic','PETsignClinic'},pathFig)
displayKaplanMeierRF_HN(pathRF,nameOutcomes,{'CTsignClinic','CTsignClinic','CTsignClinic'},pathFig)
displayKaplanMeierRF_HN(pathRF,nameOutcomes,{'CTorigsignClinic','CTorigsignClinic','CTorigsignClinic'},pathFig)
fprintf('DONE!\n'), toc

% G.3 DISPLAYING KAPLAN-MEIER CURVES FOR 3 RISK GROUPS (LOW-MEDIUM-HIGH: p < 0.33, p >= 0.33 & p < 0.66, p > 0.66)
tic, fprintf('\n--> DISPLAYING KAPLAN-MEIER CURVES FOR 3 RISK GROUPS (LOW-MEDIUM-HIGH: p < 0.33, p >= 0.33 & p < 0.66, p > 0.66) ...')
displayKaplanMeierRF_HN(pathRF,nameOutcomes,{'clinic','clinic','clinic'},pathFig,[33,66])
displayKaplanMeierRF_HN(pathRF,nameOutcomes,{'VOLclinic','VOLclinic','VOLclinic'},pathFig,[33,66])
displayKaplanMeierRF_HN(pathRF,nameOutcomes,{'PETclinic','PETclinic','PETclinic'},pathFig,[33,66])
displayKaplanMeierRF_HN(pathRF,nameOutcomes,{'CTclinic','CTclinic','CTclinic'},pathFig,[33,66])
displayKaplanMeierRF_HN(pathRF,nameOutcomes,{'PETCTclinic','PETCTclinic','PETCTclinic'},pathFig,[33,66])
displayKaplanMeierRF_HN(pathRF,nameOutcomes,{'PETsignClinic','PETsignClinic','PETsignClinic'},pathFig,[33,66])
displayKaplanMeierRF_HN(pathRF,nameOutcomes,{'CTsignClinic','CTsignClinic','CTsignClinic'},pathFig,[33,66])
displayKaplanMeierRF_HN(pathRF,nameOutcomes,{'CTorigsignClinic','CTorigsignClinic','CTorigsignClinic'},pathFig,[33,66])
fprintf('DONE!\n'), toc

time = toc(tStart);
fprintf('\n\nTOTAL TIME FOR RESULTS -- PART G: %f seconds\n',time)
fprintf('-------------------------------------------------------------------------------------')
% -------------------------------------------------------------------------



% *********************** SUPPLEMENTARY DATA  ***********************
tStart = tic;
fprintf('\n\n\n********************* SUPPLEMENTARY DATA **********************')
cd(pathRF), load('training'), load('testing'), load('permTrain'), load('permTest')
cd(pathWORK), mkdir('SUPP_DATA'), cd('SUPP_DATA'), pathSUPP = pwd; 
save('training','training'), save('testing','testing'), save('permTrain','permTrain'), save('permTest','permTest'), cd(pathWORK)

% SD.1 COMPUTING AND TESTING RANDOM FORESTS FOR RADIOMICS VARIABLES ONLY
tic, fprintf('\n--> ESTIMATNG, COMPUTING AND TESTING RANDOM FORESTS FOR RADIOMICS VARIABLES ONLY ...')
nSplit = 10; suppData_RadiomicsRF(pathSUPP,nBoot,nSplit,testSplit,testCost,[seeds(3),seeds(4)],matlabPATH,pathFig)
fprintf('DONE!\n'), toc

% SD2. CALCULATING FEATURE IMPORTANCE IN FINAL RANDOM FORESTS FOR EACH OUTCOME
tic, fprintf('\n--> CALCULATING FEATURE IMPORTANCE IN FINAL RANDOM FORESTS FOR EACH OUTCOME ...')
rng(seed), seeds = ceil(1000000000*rand(5,1));
nPerms = 100; calcAllFeatureImportanceRF_SUPP(pathRF,nPerms,seeds(5),matlabPATH)
plotFeatureImportanceRF_HN(pathRF,pathFig)
fprintf('DONE!\n'), toc

time = toc(tStart);
fprintf('\n\nTOTAL TIME FOR SUPPLEMENTARY DATA: %f seconds\n',time)
fprintf('-------------------------------------------------------------------------------------')
% -------------------------------------------------------------------------




timeAll = toc(timeStartAll);
fprintf('\n\n\n---------> TOTAL TIME FOR THIS STUDY: %f hours',timeAll/3600)
fprintf('\n**************************************************************************\n')
fprintf('------------------------------- THE END ----------------------------------\n')
fprintf('**************************************************************************\n\n')
% **************************************************************************
% ------------------------------- THE END ----------------------------------
% **************************************************************************


% TO BE DONE AFTER SUBMISSION:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TO-DO LIST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) Parallelize the RESULTS - B part the same way as with random forests:
%    1 core per outcome/fSet process
% 2) Parallelize the correction of patient names
% 3) Produce code to automatically save Figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%