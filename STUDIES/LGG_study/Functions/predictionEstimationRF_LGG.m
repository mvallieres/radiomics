function [auc,sensitivity,specificity,accuracy] = predictionEstimationRF_LGG(X,Y,cat,nTrees,seed,testCost)
% -------------------------------------------------------------------------
% function [auc,sensitivity,specificity,accuracy] = predictionEstimation_RF(X,Y,cat,nTrees,seed,testCost)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function performs prediction estimation using out-of-bag estimates
% of random forests.
% -------------------------------------------------------------------------
% INPUTS:
% 1. X: Array of size [nPatient X nFeatures] in "Table" format.
%       --> Ex: See dataStruct.mat created in organizeRF_LGG.m
% 2. Y: Column vector of size [nPatient X 1] specifying the outcome status 
%       (1 or 0) for all instances.
% 3. cat: Vecor of size [nFeatures X 1] specifying if each of the features
%         in 'X' are either continuous (0) or categorical (1).
% 4. nTrees: Number of decision-trees in the random forests.
%            --> Ex: 500
% 5. seed: Numerical number to use as seed for random-forest construction.
%          --> Ex: 54288
% 6. testCost: Vector of cost values (emphasis on instances of the minority 
%              class) to be tested in the random forest construction. The
%              one providing the highest performance on out-of-bag
%              estimations will eventually be chosen.
% -------------------------------------------------------------------------
% OUTPUTS: 
% 1. auc: Area under the receiving operating characteristic curve on
%         out-of-bag estimates.
% 2. sensitivity: Sensitivity of prediction on out-of-bag estimates.
% 3. specificity: Specificity of prediction on out-of-bag estimates.
% 4. accuracy: Accuracy of prediction on out-of-bag estimates.
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

cost = findBestCost(nTrees,X,Y,cat,testCost,seed);
rng(seed), RF = TreeBagger(nTrees,X,Y,'OOBPrediction','on','Cost',[0,1/cost;1,0],'CategoricalPredictors',cat);
[prob] = oobPredict_table(RF,X); % Probability that y = 1 for all out-of-bag estimates
[auc,sensitivity,specificity,accuracy] = calcPerformMetrics(prob,Y,0.5);

end