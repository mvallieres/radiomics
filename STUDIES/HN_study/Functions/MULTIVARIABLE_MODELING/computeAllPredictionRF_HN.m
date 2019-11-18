function computeAllPredictionRF_HN(pathModels,pathResults,training,param,maxOrder,nBoot,imbalance,seed)
% -------------------------------------------------------------------------
% function computeAllPrediction_HN(pathModels,pathResults,training,param,maxOrder,nBoot,imbalance,batchNum)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function computes prediction performance estimation for a given 
% feature set type, and for all model orders of all experiments with 
% different degrees of freedom. See ref. [1,2] for more details.
% -------------------------------------------------------------------------
% REFERENCE:
% [1] Vallieres, M. et al. (2015). FDG-PET/CT radiomics models for the 
%     early prediction of different tumour outcomes in head and neck cancer.
%     The Journal of Nuclear Medicine, aa(bb), xxx-yyy. 
%     doi:
% [2] Vallieres, M. et al. (2015). A radiomics model from joint FDG-PET and 
%     MRI texture features for the prediction of lung metastases in soft-tissue 
%     sarcomas of the extremities. Physics in Medicine and Biology, 60(14), 
%     5471-5496. doi:10.1088/0031-9155/60/14/5471
% -------------------------------------------------------------------------
% INPUTS:
% 2. pathModels: Full path to the Models folder of the corresponding experiment.
%                --> Ex: '/myProject/WORKSPACE/COHORT-BASED-RESULTS/Experiment1/MODELS'
% 2. pathResults: Full path to the Results folder of the corresponding experiment.
%                --> Ex: '/myProject/WORKSPACE/COHORT-BASED-RESULTS/Experiment1/RESULTS'
% 3. training: Structure defining all parameters for the given experiment 
%              to perform. See masterScript_HN.m for more details.
% 4. param: Cell of two strings, defining 1) The feature set type/name; and
%           2) the outcome to model
% 5. maxOrder: Integer specifying the maximal model order to construct.
%              --> Ex: 10
% 6. nBoot: Number of bootstrap samples to use.
%           --> Ex: 100
% 7. imbalance: String specifying the type of imbalance-adjustement strategy
%               employed. Either 'IABR' for imbalance-adjusted bootstrap
%               resampling (see ref.[1]), or 'IALR' for imbalance-adjusted
%               logistic regression (see ref.[2]).
%               --> Ex: 'IALR'
% 8. batchNum: (optional input). If present, integer that specifies the
%              batch number for parallelization purposes.
%              --> Ex: 6
% -------------------------------------------------------------------------
% OUTPUTS: Prediction performance results are saved in a folder named 'RESULTS' in the
% corresponding folder of each experiment (e.g. 'Experiment1', 'Experiment2', etc.)
% -------------------------------------------------------------------------
% AUTHOR(S): Martin Vallieres <mart.vallieres@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: July 2015
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
warning off

cd(pathModels)
fSetName = param{1};
outcomeName = param{2};
outcome = training.(outcomeName).outcome;
models = load(['MODELS_',fSetName,'_',outcomeName]); models = struct2cell(models); models = models{1};

nOrders = numel(fieldnames(models));
if nOrders < maxOrder
    maxOrder = nOrders;
end

tic
fprintf(['\n --> COMPUTING PREDICTION PERFORMANCE (MODEL ORDERS OF 1 to % u) FOR  "',fSetName,'" FEATURE SET AND "',outcomeName,'" OUTCOME ... '],maxOrder)
for j = 1:maxOrder
    orderName = ['Order',num2str(j)];
    data = models.(orderName).Data;
    categories = models.(orderName).Categories;
    [orderResults] = predictionPerformanceEstimationRF_HN(data,categories,outcome,nBoot,imbalance,seed);
    results.(orderName) = orderResults;
    results.(orderName).Data = models.(orderName).Data;
    results.(orderName).Categories = models.(orderName).Categories;
    results.(orderName).Name = models.(orderName).Name;
end
cd(pathResults), save(['RESULTS_',fSetName,'_',outcomeName],'results')
fprintf('DONE!\n'), toc

cd(startpath)
end
