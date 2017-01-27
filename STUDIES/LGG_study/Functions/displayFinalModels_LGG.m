function displayFinalModels_LGG(pathExperiments,nameOutcomes,fSetNames,pathFig)
% -------------------------------------------------------------------------
% function displayFinalModels_LGG(pathExperiments,nameOutcomes,fSetNames,pathFig)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function displays the final models for all analyzed outcomes in the
% MRIliver study, in terms of: multivariable model response, model variables, and 
% prediction performance estimation. See ref.[1] for more details.
% -------------------------------------------------------------------------
% REFERENCE:
% [1] Vallieres, M. et al. (2015). A radiomics model from joint FDG-PET and 
%     MRI texture features for the prediction of lung metastases in soft-tissue 
%     sarcomas of the extremities. Physics in Medicine and Biology, 60(14), 
%     5471-5496. doi:10.1088/0031-9155/60/14/5471
% -------------------------------------------------------------------------
% INPUTS:
% 1. pathExperiments: Full path to the directory containing all experiments.
%                     --> Ex: '/myProject/WORKSPACE/LOGISTIC_REGRESSION'
% 2. fSetNames: Cell of strings specifying the name of the type of feature 
%               set analyzed.
%               --> Ex: {'T1W_T2W','T1W_T2F','T1CE_T2W','T1CE_T2F'}
% 3. nameOutcomes: Cell of strings specifying the outcome names to analyze.
%                  --> Ex: {'nonIDH1','IDHcodel','progression','lowGrade'}
% 4. pathFig: (optional).  Full path to where figure is saved without
%             displaying it. Put '' for displaying the figure and not 
%             saving it to 'pathFig' (default).
%             --> Ex: ''
% -------------------------------------------------------------------------
% AUTHOR(S): Martin Vallieres <mart.vallieres@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: January 2017
%--------------------------------------------------------------------------
% STATEMENT:
% This file is part of masterScript_LiverMRI.m, a program to detect 
% desmoplastic lesions via texture analysis of MRI images.
% --> Copyright (C) 2015-2017  Martin Vallieres
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program. If not, see <http://www.gnu.org/licenses/>.
% -------------------------------------------------------------------------

if nargin < 4
    pathFig = '';
end

startpath = pwd;
nOutcomes = length(nameOutcomes);

cd(pathExperiments), load('training')

for o = 1:nOutcomes
    cd([pathExperiments,'/FINAL_MODELS/',nameOutcomes{o},'/',fSetNames{o}])
    load('finalModel'), load('coeff'), load('modelCI'), load('response') % Variables 'finalModel', 'coeff',  'modelCI', and 'response' now present in MATLAB workspace
    finalModel.outcomeData = training.(nameOutcomes{o}).outcome;
    order = size(finalModel.Data,2);
    plotSigmoidalResponse_LGG(response,finalModel.outcomeData,modelCI,nameOutcomes{o},pathFig);
    fprintf(['\n\n --> THE FINAL MULTIVARIABLE MODEL FOR "%s" OUTCOME, "%s" FSET IS:\n\n'...
             '               g(x) =               \n'],nameOutcomes{o},fSetNames{o})
    for i = 1:order
        fprintf([num2str(coeff(i)),' X ',finalModel.Name{i},'\n'])
        fprintf('                    +               \n')
    end
    fprintf(['                   ',num2str(coeff(end)),'\n'])
    fprintf('\nWITH CORRESPONDING PREDICTION PERFORMANCE ESTIMATION:\n')
    try
        fprintf(['AUC = ',num2str(roundsd(finalModel.AUC,ceil(log10(finalModel.AUC632/roundsd(finalModel.SE_AUC,1))))),' ± ',num2str(roundsd(finalModel.SE_AUC,1)),'\n'])
        fprintf(['Sensitivity = ',num2str(roundsd(finalModel.Sensitivity,ceil(log10(finalModel.Sensitivity/roundsd(finalModel.SE_Sensitivity,1))))),' ± ',num2str(roundsd(finalModel.SE_Sensitivity,1)),'\n'])
        fprintf(['Specificity = ',num2str(roundsd(finalModel.Specificity,ceil(log10(finalModel.Specificity/roundsd(finalModel.SE_Specificity,1))))),' ± ',num2str(roundsd(finalModel.SE_Specificity,1)),'\n'])
    end
    % Calculation of prediction confidence probability
    indPos = find(finalModel.outcomeData); indNeg = find(~finalModel.outcomeData);
    prob = 1./(1 + exp(-response));
    probConfidence = (sum(prob(indPos)) + sum(1 - prob(indNeg)))/numel(prob); finalModel.probConfidence = probConfidence; save('finalModel','finalModel')
    fprintf('\n*** WITH CORRESPONDING PROBABILITY OF PREDICTION CONFIDENCE: %.2f ***\n',probConfidence);
    fprintf('\n\n')
end

cd(startpath)
end