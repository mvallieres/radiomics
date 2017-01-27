function computeAllPrediction_batchLGG(pathExperiments,maxOrder,nBoot,imbalance,nBatch,matlabPATH,seed)
% -------------------------------------------------------------------------
% function computeAllPrediction_batchLGG(pathExperiments,maxOrder,nBoot,imbalance,nBatch,matlabPATH,seed)
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
% 1. pathExperiments: Full path to where all experiments need to be
%                     performed.
%                     --> Ex: /myProject/WORKSPACE/LOGISTIC_REGRESSION
% 2. maxOrder: Integer specifying the maximal model order to construct.
%              --> Ex: 10
% 3. nBoot: Number of bootstrap samples to use.
%           --> Ex: 100
% 4. imbalance: String specifying the type of imbalance-adjustement strategy
%               employed. Either 'IABR' for imbalance-adjusted bootstrap
%               resampling (see ref.[1]), or 'IALR' for imbalance-adjusted
%               logistic regression (formal reference to come).
%               --> Ex: 'IALR'
% 5. nBatch: Number of parallel batch.
%            --> Ex: 8
% 6. matlabPATH: Full path to the MATLAB excutable on the system.
%                --> 'matlab' if a symbolic link to the matlab executable
%                     was previously created.
% 7. seed: Random generator seed for reproducibility of experiments.
% -------------------------------------------------------------------------
% OUTPUTS: Prediction performance results are saved in a folder named 'RESULTS' in the
% corresponding folder of 'pathExperiments'.
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

% INITIALIZATON
time = 60; % Number of seconds to wait before checking if parallel computations are done
cd(pathExperiments), load('training')
pathModels = fullfile(pwd,'MODELS'); mkdir('RESULTS'), cd('RESULTS'), pathResults = pwd; 
mkdir('batchLog_Results'), cd('batchLog_Results'), pathBatch = pwd;
nameOutcomes = fieldnames(training); nOutcomes = numel(nameOutcomes);
for o = 1:nOutcomes
    outcomes.(nameOutcomes{o}) = training.(nameOutcomes{o}).outcome;
end
setNames = fieldnames(training.(nameOutcomes{1}).text);
[param] = batchExperiments(setNames,outcomes,nBatch); nBatch = length(param);

% PRODUCE BATCH COMPUTATIONS
save('workspace','pathModels','pathResults','training','param','maxOrder','nBoot','imbalance','seed'), pause(5);
for i = 1:nBatch
    nameScript = ['batch',num2str(i),'_script.m'];
    fid = fopen(nameScript,'w');
    fprintf(fid,'load(''workspace'')\n');
    for j = 1:numel(param{i})
        fprintf(fid,['computeAllPrediction_LGG(pathModels,pathResults,training,param{',num2str(i),'}{',num2str(j),'},maxOrder,nBoot,imbalance,seed)\n']);
    end
    fprintf(fid,['system(''touch batch',num2str(i),'_end'');\n']);
    fprintf(fid,'clear all');
    fclose(fid);
    system([matlabPATH,' -nojvm -nodisplay -nodesktop -nosplash < ',nameScript,' >& ',nameScript(1:end-1),'log &']);
end

% WAITING LOOP
waitBatch(pathBatch,time,nBatch)
delete('workspace.mat')

cd(startpath)
end