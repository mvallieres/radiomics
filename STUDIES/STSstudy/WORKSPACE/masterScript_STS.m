% *************************************************************************************
% * DESCRIPTION:                                                                      *
% * This script computes the experiments performed in ref. [1].                       *
% * --------------------------------------------------------------------------------- *
% * REFERENCE:                                                                        *
% * [1] Vallieres, M. et al. (2015). A radiomics model from joint FDG-PET and         * 
% *     MRI texture features for the prediction of lung metastases in soft-tissue     * 
% *     sarcomas of the extremities. Physics in Medicine and Biology, 60(14),         * 
% *     5471-5496. doi:10.1088/0031-9155/60/14/5471                                   * 
% * --------------------------------------------------------------------------------- *
% * AUTHOR(S): Martin Vallieres <mart.vallieres@gmail.com>                            *
% * --------------------------------------------------------------------------------- *
% * HISTORY:                                                                          *
% * - Creation: May 2015                                                              *
% * --------------------------------------------------------------------------------- *
% * STATEMENT:                                                                        *
% * This file is part of <https://github.com/mvallieres/radiomics/>,                  *
% * a package providing MATLAB programming tools for radiomics analysis.              *
% * --> Copyright (C) 2015  Martin Vallieres                                          *
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

fprintf('\n')
help masterScript_STS
fprintf('\n')
warning off

% INITIALIZATION
pathWORK = pwd;
roiNumb = load('contour_Mass'); roiNumb = struct2cell(roiNumb); roiNumb = roiNumb{1};
outcome = load('outcome'); outcome = struct2cell(outcome); outcome = outcome{1};
nPatient = numel(outcome);

% TEXTURE EXTRACTION PARAMETERS AND DEGREES OF FREEDOM
MRIinv_cell = {'NoInv','Inv'};
MRIweight_mat = [1/4,1/3,1/2,2/3,3/4];
R_mat = [1/2,2/3,1,3/2,2];
scale_cell = {'pixelW',1,2,3,4,5};
algo_cell = {'Equal','Lloyd'};
Ng_mat = [8,16,32,64];
paramSEP = {R_mat,scale_cell,algo_cell,Ng_mat};                           % NOTE: paramSEP and paramFUS must always be of the same format and size for SEPARATE 
paramFUS = {MRIinv_cell,MRIweight_mat,R_mat,scale_cell,algo_cell,Ng_mat}; % and FUSED SCANS as presented here, with the same ordering of different extraction parameters.
baselineSEP = [3,1,2,3]; % As defined in ref. [1]
baselineFUS = [1,3,3,1,2,3]; % As defined in ref. [1]
freedomSEP = [unique(perms([0 0 0 0]),'rows');unique(perms([0 0 0 1]),'rows');...
              unique(perms([0 0 1 1]),'rows');unique(perms([0 1 1 1]),'rows');...
              unique(perms([1 1 1 1]),'rows')];
freedomFUS = [unique(perms([0 0 0 0 0 0]),'rows');unique(perms([0 0 0 0 0 1]),'rows');...
              unique(perms([0 0 0 0 1 1]),'rows');unique(perms([0 0 0 1 1 1]),'rows');...
              unique(perms([0 0 1 1 1 1]),'rows');unique(perms([0 1 1 1 1 1]),'rows');...
              unique(perms([1 1 1 1 1 1]),'rows')];

% MULTIVARIABLE ANALYSIS PARAMETERS
nBoot = 1000;
alpha = 0.5; delta = 0.5;
setSize = 25;
fSetName = {'PET','SEPARATE','FUSED'};
maxOrder = 10;
          


% ************** START COMPUTATION **************
          
% 1. READ DATA DOWNLOADED FROM THE TCIA WEBSITE (http://dx.doi.org/10.7937/K9/TCIA.2015.7GO2GSKS)
fprintf('\n\n*********************** ORGANIZING AND PROCESSING DICOM DATA FROM TCIA WEBSITE ***********************\n')
readAllDICOM_STS([pathWORK,'/Soft-tissue-Sarcoma'],nPatient)
 
% 2. COMPUTE NON-TEXTURE FEATURES
fprintf('\n\n*********************** COMPUTING NON-TEXTURE FEATURES ***********************\n')
calcAllNonTextureFeatures_STS(pathWORK,nPatient,roiNumb,outcome)

% 3. COMPUTE AND ORGANIZING ALL TEXTURE FEATURES
mkdir('TEXTURES'), fprintf('\n')
calcAllSeparateTextures_STS(pathWORK,nPatient,roiNumb,R_mat,scale_cell,algo_cell,Ng_mat)
organizeSeparateTextures_STS(pathWORK,nPatient,R_mat,scale_cell,algo_cell,Ng_mat)
calcAllFusedTextures_STS(pathWORK,nPatient,roiNumb,MRIinv_cell,MRIweight_mat,R_mat,scale_cell,algo_cell,Ng_mat)
organizeFusedTextures_STS(pathWORK,nPatient,MRIinv_cell,MRIweight_mat,R_mat,scale_cell,algo_cell,Ng_mat)



% 5. PERFORM FEATURE SET REDUCTION
mkdir('FSET')
nonTextures = load('nonTextures'); nonTextures = struct2cell(nonTextures); nonTextures = nonTextures{1};
PET  = load('textures_PET');  PET  = struct2cell(PET);  PET  = PET{1};
T1   = load('textures_T1');   T1   = struct2cell(T1);   T1   = T1{1};
T2FS = load('textures_T2FS'); T2FS = struct2cell(T2FS); T2FS = T2FS{1};
PET_T1   = load('textures_PET_T1');   PET_T1   = struct2cell(PET_T1);   PET_T1   = PET_T1{1};
PET_T2FS = load('textures_PET_T2FS'); PET_T2FS = struct2cell(PET_T2FS); PET_T2FS = PET_T2FS{1};

tic
fprintf('\n\nFinding the path to the ''MINE.jar'' application on the system ... ')
pathMINE = findMINE('Linux');
fprintf('DONE\n')
toc

% For Feature set 1: PET
fprintf('\n*********************** PERFORMING FEATURE SET REDUCTION FOR ''PET'' FEATURE SET ***********************\n')
calcAllFeatureSets_STS(pathWORK,pathMINE,fSetName{1},outcome,setSize,nonTextures,{PET},{'PET'},paramSEP,freedomSEP,baselineSEP,alpha,delta,nBoot)

% For Feature set 2: PET, T1, T2FS
fprintf('\n*********************** PERFORMING FEATURE SET REDUCTION FOR ''SEPARATE'' FEATURE SET **********************\n')
calcAllFeatureSets_STS(pathWORK,pathMINE,fSetName{2},outcome,setSize,nonTextures,{PET,T1,T2FS},{'PET','T1','T2FS'},paramSEP,freedomSEP,baselineSEP,alpha,delta,nBoot)

% For Feature set 3: PET_T1, PET_T2FS
fprintf('\n*********************** PERFORMING FEATURE SET REDUCTION FOR ''FUSED'' FEATURE SET ***********************\n')
calcAllFeatureSets_STS(pathWORK,pathMINE,fSetName{3},outcome,setSize,nonTextures,{PET_T1,PET_T2FS},{'PET_T1','PET_T2FS'},paramFUS,freedomFUS,baselineFUS,alpha,delta,nBoot)



% 6. PERFORM FEATURE SET SELECTION
mkdir('MODELS')
fprintf('\n\nFinding the path to''fastAUC.cpp'' on the system --> COMPILATION ... ')
try 
    compileFastAUC('Linux')
    fprintf('DONE\n') 
catch
    fprintf('FAILED (AUC computations will be slower)\n')
end
 
% For Feature set 1: PET
fprintf('\n*********************** PERFORMING FEATURE SELECTION FOR ''PET'' FEATURE SET ***********************\n')
computeAllModelChoice_STS(pathWORK,fSetName{1},outcome,freedomSEP,maxOrder,nBoot)

% For Feature set 2: PET, T1, T2FS
fprintf('\n*********************** PERFORMING FEATURE SELECTION  FOR ''SEPARATE'' FEATURE SET ***********************\n')
computeAllModelChoice_STS(pathWORK,fSetName{2},outcome,freedomSEP,maxOrder,nBoot)

% For Feature set 3: PET_T1, PET_T2FS
fprintf('\n*********************** PERFORMING FEATURE SELECTION  FOR ''FUSED'' FEATURE SET ***********************\n')
computeAllModelChoice_STS(pathWORK,fSetName{3},outcome,freedomFUS,maxOrder,nBoot)



% 7. PERFORM PREDICTION PERFORMANCE ESTIMATION
mkdir('RESULTS'), fprintf('\n')

% For Feature set 1: PET
fprintf('\n*********************** PERFORMING PREDICTION PERFORMANCE ESTIMATION FOR ''PET'' FEATURE SET ***********************\n')
computeAllPrediction_STS(pathWORK,fSetName{1},outcome,freedomSEP,maxOrder,nBoot)

% For Feature set 2: PET, T1, T2FS
fprintf('\n*********************** PERFORMING PREDICTION PERFORMANCE ESTIMATION FOR ''SEPARATE'' FEATURE SET ***********************\n')
computeAllPrediction_STS(pathWORK,fSetName{2},outcome,freedomSEP,maxOrder,nBoot)

% For Feature set 3: PET_T1, PET_T2FS
fprintf('\n*********************** PERFORMING PREDICTION PERFORMANCE ESTIMATION FOR ''FUSED'' FEATURE SET ***********************\n')
computeAllPrediction_STS(pathWORK,fSetName{3},outcome,freedomFUS,maxOrder,nBoot)

% Finding the best combinations of model order and texture extraction parameter degree of freedom for all feature set types
groupExperiments_STS(pathWORK,fSetName{1},freedomSEP,maxOrder)
groupExperiments_STS(pathWORK,fSetName{2},freedomSEP,maxOrder)
groupExperiments_STS(pathWORK,fSetName{3},freedomFUS,maxOrder)



% 8. CHOICE OF BEST PARSIMONIOUS MODEL (requires user input)
fprintf('\n')
plotPredictionResults_STS([pathWORK,'/RESULTS'],fSetName,{'AUC632','Sensitivity632','Specificity632'},maxOrder)
while 1
    set = input(['\nWhich feature set provides the best parsimonious model? \n' ...
                 '--> For the ',fSetName{1},' feature set, type ''1'' and press ENTER \n' ...
                 '--> For the ',fSetName{2},' feature set, type ''2'' and press ENTER \n' ...
                 '--> For the ',fSetName{3},' feature set, type ''3'' and press ENTER \n' ...
                 'ANSWER: ']);
    fprintf('\n')
    if isnumeric(set) && (set == 1 || set == 2 || set == 3)
        break
    end
end
while 1
    order = input(['Which model order of the ',fSetName{set},' feature set provides the best parsimonious model? \n' ...
                   '--> Type a number between 1 to ',num2str(maxOrder),' and press ENTER \n' ...
                   'ANSWER: ']);
    fprintf('\n')
    if isnumeric(order) && order <= 10 && order >= 1
        break
    end
end
cd([pathWORK,'/RESULTS'])
results = load(['RESULTS_',fSetName{set},'_BEST']); results = struct2cell(results); results = results{1}; 
finalModel = results.(['Order',num2str(order)]); 
cd(pathWORK), mkdir('FINAL_MODEL'), cd('FINAL_MODEL'), save('finalModel', 'finalModel')



% 9. COMPUTING THE LOGISTIC REGRESSION COEFFICIENTS AND BOOTSTRAP CONFIDENCE INTERVALS OF THE FINAL MODEL
fprintf('\n\nCOMPUTING THE LOGISTIC REGRESSION COEFFICIENTS OF THE FINAL MODEL ... ')
tic
[coeff,response,modelCI] = computeModelCoefficients(finalModel.Data,outcome,'IABR');
save('coeff','coeff'), save('response','response'), save('modelCI','modelCI')
fprintf('DONE\n')
toc



% 10. DISPLAYING THE FINAL MODEL AND CORRESPONDING PREDICTION PERFORMANCE ESTIMATION
plotSigmoidalResponse(response,outcome,modelCI,'LungMets')
fprintf(['\n\n\n --> THE FINAL MULTIVARIABLE MODEL IS:\n\n'...
         '               g(x) =               \n'])
for i = 1:order
    fprintf([num2str(coeff(i)),' X ',finalModel.Name{i},'\n'])
    fprintf('                    +               \n')
end
fprintf(['                   ',num2str(coeff(end)),'\n'])
fprintf('\nWITH CORRESPONDING PREDICTION PERFORMANCE ESTIMATION:\n')
fprintf(['AUC = ',num2str(roundsd(finalModel.AUC632,ceil(log10(finalModel.AUC632/roundsd(finalModel.SE_AUC632,1))))),' ± ',num2str(roundsd(finalModel.SE_AUC632,1)),'\n'])
fprintf(['Sensitivity = ',num2str(roundsd(finalModel.Sensitivity632,ceil(log10(finalModel.Sensitivity632/roundsd(finalModel.SE_Sensitivity632,1))))),' ± ',num2str(roundsd(finalModel.SE_Sensitivity632,1)),'\n'])
fprintf(['Specificity = ',num2str(roundsd(finalModel.Specificity632,ceil(log10(finalModel.Specificity632/roundsd(finalModel.SE_Specificity632,1))))),' ± ',num2str(roundsd(finalModel.SE_Specificity632,1)),'\n'])
fprintf('\n')
cd(pathWORK)
