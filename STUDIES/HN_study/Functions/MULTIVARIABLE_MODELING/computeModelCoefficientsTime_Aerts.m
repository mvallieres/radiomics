function computeModelCoefficientsTime_Aerts(pathExperimentsTime,seed)
% -------------------------------------------------------------------------
% function computeModelCoefficients_batchHN(pathExperiments,nExp,fSetNames,imbalance,matlabPATH)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function computes the final logistic regression coefficients and
% bootstrap confidence intervals of the final models obtained for all
% outcomes analyzed in the HN study. See ref.[1] for more details.
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
% 3. imbalance: String specifying the type of imbalance-adjustement strategy
%               employed. Either 'IABR' for imbalance-adjusted bootstrap
%               resampling (see ref.[1]), or 'IALR' for imbalance-adjusted
%               logistic regression (see ref.[2]).
%               --> Ex: 'IALR'
% 4. nBatch: Number of parallel batch.
%            --> Ex: 8
% 5. matlabPATH: Full path to the MATLAB excutable on the system.
%                --> 'matlab' if a symbolic link to the matlab executable
%                     was previously created.
% -------------------------------------------------------------------------
% OUTPUTS: Final coefficients, model response and confidence intervals 
%          saved in the corresponding ../experiment/outcome/fSetName
%          folder.
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


% INITIALIZATIONS
cd(pathExperimentsTime), load('training')
cd('FINAL_MODELS'), mkdir('DeathSign'), cd('DeathSign'), mkdir('CT'), cd('CT'), pathFinalModel = pwd;
finalModel.Data = training.Death.sign.CT;
finalModel.Name = {'CT_Energy';'CT_Compactness';'CT_GLN';'CT_GLN_HLH'};
finalModel.Order = 4;
finalModel.outcome = 'Death';
save('finalModel','finalModel')

% COMPUTATIONS
[coeff,response,modelCI,medianHR] = computeModelCoefficientsTime_HN(finalModel.Data,training.Death.timeToEvent,1 - training.Death.outcome,seed);
save('coeff','coeff'), save('response','response'), save('modelCI','modelCI'), save('medianHR','medianHR')

cd(startpath)
end