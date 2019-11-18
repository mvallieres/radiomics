function autoChoiceBestModels_HN(pathExperimentsBinary,pathExperimentsTime,fSetNames,nameOutcomes,metric,maxOrder,pathFig)
% -------------------------------------------------------------------------
% function chooseBestModels_HN(pathWORK,nExp,fSetName,nameOutcomes,maxOrder)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function requires user input to choose the set types and model orders
% providing the best parsimonious models for all outcomes analyzed in the
% HN study. See ref.[1] for more details.
% -------------------------------------------------------------------------
% REFERENCE:
% [1] Vallieres, M. et al. (2015). FDG-PET/CT radiomics models for the 
%     early prediction of different tumour outcomes in head and neck cancer.
%     The Journal of Nuclear Medicine, aa(bb), xxx-yyy. 
%     doi:
% -------------------------------------------------------------------------
% INPUTS:
% 1. pathExperiments: Full path to the directory containing all experiments.
%                     --> Ex: '/myProject/WORKSPACE/CV-BASED_RESULTS'
% 2. nExp: Numerical value specifying the number of experiments to analyze.
%          --> Ex: 10
% 3. fSetName: Cell of strings specifying the name of the type of feature 
%              set analyzed.
%              --> Ex: {'PET','CT','SEPARATE','FUSED'}
% 4. nameOutcomes: Cell of strings specifying the outcome names to analyze.
%                  --> Ex: {'Failure','Locoregional','Distant','Death'}
% 5. maxOrder: Integer specifying the maximal multivariable model order.
%              --> Ex: 10
% -------------------------------------------------------------------------
% OUTPUTS: Final prediction models saved in a folder named 'FINAL_MODELS/nameOutcome/fSetName'
%          in the given experiment folder.
% -------------------------------------------------------------------------
% AUTHOR(S): Martin Vallieres <mart.vallieres@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: March 2016
%--------------------------------------------------------------------------
% STATEMENT:
% This file is part of <https://github.com/mvallieres/radiomics/>, 
% a package providing MATLAB programming tools for radiomics analysis.
% --> Copyright (C) 2015  Martin Vallieres
%
%    This package is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This package is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this package.  If not, see <http://www.gnu.org/licenses/>.
% -------------------------------------------------------------------------

startpath = pwd;
cd(pathExperimentsBinary)
nOutcomes = length(nameOutcomes);
nType = length(fSetNames);

n = 0;
%order = ones(1,nOutcomes*nType)*3; % All order 3
order = [8,3,3,6,3,3,4,3,6];

cd(pathExperimentsTime), pathExperimentTime = pathExperimentsTime;
cd(pathExperimentTime), mkdir('FINAL_MODELS'), cd('FINAL_MODELS'), pathFinalModelsTime = pwd;
cd(pathExperimentsBinary), pathExperimentBinary = pathExperimentsBinary;
cd(pathExperimentBinary), mkdir('FINAL_MODELS'), cd('FINAL_MODELS'), pathFinalModelsBinary = pwd;
for o = 1:nOutcomes
    cd(pathFinalModelsTime), mkdir(nameOutcomes{o}), cd(fullfile(pathFinalModelsTime,nameOutcomes{o})), pathOutcomeTime = pwd;
    cd(pathFinalModelsBinary), mkdir(nameOutcomes{o}), cd(fullfile(pathFinalModelsBinary,nameOutcomes{o})), pathOutcomeBinary = pwd;
    plotPredictionResults_HN(fullfile(pathExperimentBinary,'RESULTS'),nameOutcomes{o},fSetNames,metric,maxOrder,pathFig)
    for f = 1:nType
        n = n + 1;
        cd(pathOutcomeTime), mkdir(fSetNames{f}), cd(fullfile(pathOutcomeTime,fSetNames{f})), pathFinalTime = pwd;
        cd(pathOutcomeBinary), mkdir(fSetNames{f}), cd(fullfile(pathOutcomeBinary,fSetNames{f})), pathFinalBinary = pwd;
        cd(fullfile(pathExperimentBinary,'RESULTS'))
        results = load(['RESULTS_',fSetNames{f},'_',nameOutcomes{o}]); results = struct2cell(results); results = results{1};
        finalModel = results.(['Order',num2str(order(n))]);
        finalModel.Order = order(n);
        finalModel.outcome = nameOutcomes{o};
        cd(pathFinalBinary), save('finalModel','finalModel')
        cd(pathFinalTime), save('finalModel','finalModel')
    end
end

cd(startpath)
end