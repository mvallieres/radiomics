function computeAllModelChoice_LGG(pathSet,pathModels,training,param,maxOrder,nBoot,imbalance,seed)
% -------------------------------------------------------------------------
% function computeAllModelChoice_LGG(pathSet,pathModels,training,param,maxOrder,nBoot,imbalance,seed)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function computes feature set selection for a given feature set type 
% and for all experiments with different degrees of freedom. See ref. [1] 
% for more details.
% -------------------------------------------------------------------------
% REFERENCES:
% [1] Vallieres, M. et al. (2015). A radiomics model from joint FDG-PET and 
%     MRI texture features for the prediction of lung metastases in soft-tissue 
%     sarcomas of the extremities. Physics in Medicine and Biology, 60(14), 
%     5471-5496. doi:10.1088/0031-9155/60/14/5471
% -------------------------------------------------------------------------
% INPUTS:
% 1. pathSet: Full path to the Feature Set folder of the corresponding experiment.
%             --> Ex: '/myProject/WORKSPACE/LOGISTIC_REGRESSION/FSET'
% 2. pathModels: Full path to the Models folder of the corresponding experiment.
%                --> Ex: '/myProject/WORKSPACE/LOGISTIC_REGRESSION/MODELS'
% 3. training: Structure defining all parameters for the given experiment 
%              to perform. See organizeRadiomicsExp_LGG.m for more details.
%              --> IMPORTANT VARIABLE DEFINING ALL PARAMETERS IN TRAINING
% 4. param: Cell of two strings, defining 1) The feature set type/name; and
%           2) the outcome to model in the 'training' structure.
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
% OUTPUTS: Mutivariable models are saved in 'pathMODELS'.
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

cd(pathSet)
fSetName = param{1};
outcomeName = param{2};
outcome = training.(outcomeName).outcome;
fSet = load(['FSET_',fSetName,'_',outcomeName]); fSet = struct2cell(fSet); fSet = fSet{1};

tic
fprintf(['\n--> SELECTING FEATURES (MODEL ORDERS OF 1 to % u) FOR "',fSetName,'" FEATURE SET AND "',outcomeName,'" OUTCOME ... '],maxOrder)
[models] = featureSelection(fSet.Data,outcome,maxOrder,nBoot,fSet.Info,imbalance,seed);
cd(pathModels), save(['MODELS_',fSetName,'_',outcomeName],'models')
fprintf('DONE!\n'), toc

cd(startpath)
end
