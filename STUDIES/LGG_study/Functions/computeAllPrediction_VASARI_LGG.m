function computeAllPrediction_VASARI_LGG(pathModels,pathResults,outcomes,param,maxOrder,nBoot,imbalance,seed)
% -------------------------------------------------------------------------
% function computeAllPrediction_VASARI_LGG(pathModels,pathResults,outcomes,param,maxOrder,nBoot,imbalance,seed)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function computes prediction performance estimation for a given 
% feature set type, and for all model orders of all experiments with 
% different degrees of freedom. See ref. [1] for more details.
% -------------------------------------------------------------------------
% REFERENCE:
% [1] Vallieres, M. et al. (2015). A radiomics model from joint FDG-PET and 
%     MRI texture features for the prediction of lung metastases in soft-tissue 
%     sarcomas of the extremities. Physics in Medicine and Biology, 60(14), 
%     5471-5496. doi:10.1088/0031-9155/60/14/5471
% -------------------------------------------------------------------------
% INPUTS:
% 1. pathModels: Full path to the Models folder of the corresponding experiment.
%                --> Ex: '/myProject/WORKSPACE/VASARI/MODELS'
% 2. pathResults: Full path to the Results folder of the corresponding experiment.
%                --> Ex: '/myProject/WORKSPACE/VASARI/RESULTS'
% 3. outcomes: Structure specifying the status (1 or 0) for different
%              outcomes in LGG cancer. Contains: outcomes.nonIDH1, 
%              outcomes.IDHcodel, outcomes.progression and
%              outcomes.lowGrade.
% 4. param: Cell of two strings, defining 1) The feature set type/name; and
%           2) the outcome to model
% 5. maxOrder: Integer specifying the maximal model order to construct.
%              --> Ex: 10
% 6. nBoot: Number of bootstrap samples to use.
%           --> Ex: 100
% 7. imbalance: String specifying the type of imbalance-adjustement strategy
%               employed. Either 'IABR' for imbalance-adjusted bootstrap
%               resampling (see ref.[1]), or 'IALR' for imbalance-adjusted
%               logistic regression (formal reference to come).
%               --> Ex: 'IALR'
% 8. seed: Numerical number to use as seed for bootstrapping experiment
%          --> Ex: 54288
% -------------------------------------------------------------------------
% OUTPUTS: Mutivariable models are saved in a folder named 'RESULTS'.
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
warning off

cd(pathModels)
fSetName = param{1};
outcomeName = param{2};
outcome = outcomes.(outcomeName);
models = load(['MODELS_',fSetName,'_',outcomeName]); models = struct2cell(models); models = models{1};

tic
fprintf(['\n --> COMPUTING PREDICTION PERFORMANCE (MODEL ORDERS OF 1 to % u) FOR  "',fSetName,'" FEATURE SET AND "',outcomeName,'" OUTCOME ... '],maxOrder)
for j = 1:maxOrder
    orderName = ['Order',num2str(j)];
    data = models.(orderName).Data;
    [orderResults] = predictionPerformanceEstimation(data,outcome,nBoot,imbalance,seed);
    results.(orderName) = orderResults;
    results.(orderName).Data = models.(orderName).Data;
    results.(orderName).Name = models.(orderName).Name;
end
cd(pathResults), save(['RESULTS_',fSetName,'_',outcomeName],'results')
fprintf('DONE!\n'), toc

cd(startpath)
end
