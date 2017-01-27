function computeAllModelChoice_batchVASARI_LGG(pathExperiments,maxOrder,nBoot,imbalance,nBatch,matlabPATH,seed)
% -------------------------------------------------------------------------
% function computeAllModelChoice_batchVASARI_LGG(pathExperiments,maxOrder,nBoot,imbalance,nBatch,matlabPATH,seed)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function computes feature set selection for all feature set types 
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
% 1. pathExperiments: Full path to where all experiments need to be
%                     performed.
%                     --> Ex: /myProject/WORKSPACE/VASARI
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
% 6. matlabPATH: Full path to the MATLAB executable on the system.
%                --> 'matlab' if a symbolic link to the matlab executable
%                     was previously created.
% 7. seed: Numerical number to use as seed for bootstrapping experiment
%          --> Ex: 54288
% -------------------------------------------------------------------------
% OUTPUTS: Mutivariable models are saved in a folder named 'MODELS' in the
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
cd(pathExperiments), load('outcomes')
pathSet = fullfile(pwd,'FSET'); mkdir('MODELS'), cd('MODELS'), pathModels = pwd; 
mkdir('batchLog_Models'), cd('batchLog_Models'), pathBatch = pwd;
setNames = {'VASARI'};
[param] = batchExperiments(setNames,outcomes,nBatch); nBatch = length(param);

% PRODUCE BATCH COMPUTATIONS
save('workspace','pathSet','pathModels','outcomes','param','maxOrder','nBoot','imbalance','seed'), pause(5);
for i = 1:nBatch
    nameScript = ['batch',num2str(i),'_script.m'];
    fid = fopen(nameScript,'w');
    fprintf(fid,'load(''workspace'')\n');
    for j = 1:numel(param{i})
        fprintf(fid,['computeAllModelChoice_VASARI_LGG(pathSet,pathModels,outcomes,param{',num2str(i),'}{',num2str(j),'},maxOrder,nBoot,imbalance,seed)\n']);
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