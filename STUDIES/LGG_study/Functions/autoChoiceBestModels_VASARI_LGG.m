function autoChoiceBestModels_VASARI_LGG(pathExperiments,fSetNames,nameOutcomes,metric,maxOrder,pathFig)
% -------------------------------------------------------------------------
% function autoChoiceBestModels_VASARI_LGG(pathExperiments,fSetNames,nameOutcomes,metric,maxOrder,pathFig)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function automatically chooses (hard-coded in the function) the set 
% types and models orders providing the best parsimonious models for all
% outcoms analyzed in the LGG study.
% -------------------------------------------------------------------------
% INPUTS:
% 1. pathExperiments: Full path to the directory containing all experiments.
%                     --> Ex: '/myProject/WORKSPACE/VASARI'
% 2. fSetNames: Cell of strings specifying the name of the type of feature 
%               set analyzed.
%               --> Ex: {'VASARI'}
% 3. nameOutcomes: Cell of strings specifying the outcome names to analyze.
%                  --> Ex: {'nonIDH1','IDHcodel','progression','lowGrade'}
% 4. metric: String specifying the metric to display.
%            --> 'AUC632'
% 5. maxOrder: Integer specifying the maximal multivariable model order.
%              --> Ex: 10
% 6. pathFig: (optional).  Full path to where figure is saved without
%             displaying it. Put '' for displaying the figure and not 
%             saving it to 'pathFig' (default).
%             --> Ex: ''
% -------------------------------------------------------------------------
% OUTPUTS: Final prediction models saved in a folder named 
%          '/myProject/WORKSPACE/VASARI/FINAL_MODELS/nameOutcome/fSetName'
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
cd(pathExperiments)
nOutcomes = length(nameOutcomes);
nType = length(fSetNames);

n = 0;
order = [2,4,5,3]; % For 'nonIDH1', IDHcodel, 'progression', 'lowGrade'

cd(pathExperiments), pathExperimentBinary = pathExperiments;
cd(pathExperimentBinary), mkdir('FINAL_MODELS'), cd('FINAL_MODELS'), pathFinalModelsBinary = pwd;
for o = 1:nOutcomes
    cd(pathFinalModelsBinary), mkdir(nameOutcomes{o}), cd(fullfile(pathFinalModelsBinary,nameOutcomes{o})), pathOutcomeBinary = pwd;
    plotPredictionResults_VASARI(fullfile(pathExperimentBinary,'RESULTS'),nameOutcomes{o},fSetNames,metric,maxOrder,pathFig)
    for f = 1:nType
        n = n + 1;
        cd(pathOutcomeBinary), mkdir(fSetNames{f}), cd(fullfile(pathOutcomeBinary,fSetNames{f})), pathFinalBinary = pwd;
        cd(fullfile(pathExperimentBinary,'RESULTS'))
        results = load(['RESULTS_',fSetNames{f},'_',nameOutcomes{o}]); results = struct2cell(results); results = results{1};
        finalModel = results.(['Order',num2str(order(n))]);
        finalModel.Order = order(n);
        finalModel.outcome = nameOutcomes{o};
        cd(pathFinalBinary), save('finalModel','finalModel')
    end
end

cd(startpath)
end