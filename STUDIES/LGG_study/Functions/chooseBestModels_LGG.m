function chooseBestModels_LGG(pathExperiments,fSetNames,nameOutcomes,metric,maxOrder)
% -------------------------------------------------------------------------
% function chooseBestModels_LGG(pathExperiments,fSetNames,nameOutcomes,metric,maxOrder)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function requires user input to choose the set types and model orders
% providing the best parsimonious models for all outcomes analyzed in the
% LGG study.
% -------------------------------------------------------------------------
% INPUTS:
% 1. pathExperiments: Full path to the directory containing all experiments.
%                     --> Ex: '/myProject/WORKSPACE/LOGISTIC_REGRESSION'
% 2. fSetNames: Cell of strings specifying the name of the type of feature 
%               set analyzed.
%               --> Ex: {'T1W_T2W','T1W_T2F','T1CE_T2W','T1CE_T2F'}
% 3. nameOutcomes: Cell of strings specifying the outcome names to analyze.
%                  --> Ex: {'nonIDH1','IDHcodel','progression','lowGrade'}
% 4. metric: String specifying the metric to display.
%            --> 'AUC632'
% 5. maxOrder: Integer specifying the maximal multivariable model order.
%              --> Ex: 10
% -------------------------------------------------------------------------
% OUTPUTS: Final prediction models saved in a folder named 
%          '/myProject/WORKSPACE/LOGISTIC_REGRESSION/FINAL_MODELS/nameOutcome/fSetName'
% -------------------------------------------------------------------------
% AUTHOR(S): Martin Vallieres <mart.vallieres@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: January 2017
%--------------------------------------------------------------------------
% STATEMENT:
% This file is part of <https://github.com/mvallieres/radiomics/>, 
% a package providing MATLAB programming tools for radiomics analysis.
% --> Copyright (C) 2015-2017  Martin Vallieres
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
nOutcomes = length(nameOutcomes);
nType = length(fSetNames);

fprintf('\n')
cd(pathExperiments), mkdir('FINAL_MODELS'), cd('FINAL_MODELS'), pathFinalModels = pwd;
for o = 1:nOutcomes
    cd(pathFinalModels), mkdir(nameOutcomes{o}), cd(fullfile(pathFinalModels,nameOutcomes{o})), pathOutcome = pwd;
    plotPredictionResults_LGG(fullfile(pathExperiments,'RESULTS'),nameOutcomes{o},fSetNames,metric,maxOrder)
    fprintf('\n======== DISPLAYING PREDICTION RESULTS FOR "%s" OUTCOME, "%s" METRIC  ========\n',nameOutcomes{o},metric)
    for f = 1:nType
        cd(pathOutcome), mkdir(fSetNames{f}), cd(fullfile(pathOutcome,fSetNames{f})), pathFinal = pwd;
        while 1
            order = input(['Which model order of the ',fSetNames{f},' feature set provides the best parsimonious model? \n' ...
                           '--> Type a number between 1 to ',num2str(maxOrder),' and press ENTER \n' ...
                           'ANSWER: ']);
            fprintf('\n')
            if isnumeric(order) && order <= maxOrder && order >= 1
                break
            end
        end
        cd(fullfile(pathExperiments,'RESULTS'))
        results = load(['RESULTS_',fSetNames{f},'_',nameOutcomes{o}]); results = struct2cell(results); results = results{1};
        finalModel = results.(['Order',num2str(order)]);
        finalModel.Order = order;
        finalModel.outcome = nameOutcomes{o};
        cd(pathFinal), save('finalModel','finalModel')
    end
    close all
end

cd(startpath)
end