function cost = findBestCost(nTrees,X,Y,cat,testCost,seed)
% -------------------------------------------------------------------------
% function cost = findBestCost(nTrees,X,Y,cat,testCost,seed)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function finds the best cost (emphasis on instances of the minority 
% class) for best performance for a given random forest experiment.
% -------------------------------------------------------------------------
% INPUTS:
% 1. nTrees: Number of decision-trees in the random forests.
%            --> Ex: 500
% 2. X: Array of size [nPatient X nFeatures] in "Table" format.
%       --> Ex: See dataStruct.mat created in organizeRF_LGG.m
% 3. Y: Column vector of size [nPatient X 1] specifying the outcome status 
%       (1 or 0) for all instances.
% 4. cat: Vector of size [nFeatures X 1] specifying if each of the features
%         in 'X' are either continuous (0) or categorical (1).
% 5. testCost: Vector of cost values (emphasis on instances of the minority 
%              class) to be tested in the random forest construction. The
%              one providing the highest performance on out-of-bag
%              estimations will eventually be chosen.
% 6. seed: Numerical number to use as seed for random-forest construction.
%          --> Ex: 54288
% -------------------------------------------------------------------------
% OUTPUTS: 
% 1. cost: Numerical value
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

nCost = numel(testCost);
maxMetric = 0; cost = 1; % Default
for c = 1:nCost
    rng(seed), RF = TreeBagger(nTrees,X,Y,'OOBPrediction','on','Cost',[0,1/testCost(c);1,0],'CategoricalPredictors',cat);
    [prob] = oobPredict_table(RF,X); % Probability that y = 1 for all out-of-bag estimates
    [AUC,sensitivity,specificity,~] = calcPerformMetrics(prob,Y,0.5);
    metric = 0.5*AUC + 0.5*(1 - abs(sensitivity - specificity)); % Metric to maximize.
    if metric > maxMetric
        maxMetric = metric;
        cost = testCost(c);
    end
end

end
