function calcAllFeatureSets_batchLGG(pathMINE,pathExperiments,setSize,alpha,delta,nBoot,nBatch,matlabPATH,seed)
% -------------------------------------------------------------------------
% function calcAllFeatureSets_batchLGG(pathMINE,pathExperiments,nExperiments,setSize,alpha,delta,nBoot,nBatch,matlabPATH)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function computes feature set reduction for all feature set types
% and outcomes and for all experiments with different degrees of freedom 
% given as inputs. See ref. [1] for more details.
% -------------------------------------------------------------------------
% REFERENCES:
% [1] Vallieres, M. et al. (2015). A radiomics model from joint FDG-PET and 
%     MRI texture features for the prediction of lung metastases in soft-tissue 
%     sarcomas of the extremities. Physics in Medicine and Biology, 60(14), 
%     5471-5496. doi:10.1088/0031-9155/60/14/5471
% -------------------------------------------------------------------------
% INPUTS:
% 1. pathMINE: Full path to the MINE.jar executable on the system
%              --> Ex: '/myProject/radiomics/MultivariableModeling/MINE'
% 2. pathExperiments: Full path to where all experiments need to be
%                     performed.
%                     --> Ex: '/myProject/WORKSPACE/LOGISTIC_REGRESSION'
% 3. setSize: Size of the output feature set (typically set to 25 in ref. [1]).
%             --> Ex: 25
% 4. alpha: Numerical values specifying the coefficient of the first part of
%           the Gain equation, as defined in ref. [1].
%           --> Ex: 0.5
% 5. delta: Numerical values specifying the coefficient of the second part 
%           of the Gain equation, as defined in ref. [1] (third part is set
%           to 0 in this function).
%           --> Ex: 0.5
% 6. nBoot: Number of bootstrap samples to use.
%           --> Ex: 100
% 7. nBatch: Number of parallel batch.
%            --> Ex: 8
% 8. matlabPATH: Full path to the MATLAB executable on the system.
%                --> 'matlab' if a symbolic link to the matlab executable
%                     was previously created.
% 9. seed: Numerical number to use as seed for bootstrapping experiment
%          --> Ex: 54288
%
% See <https://github.com/mvallieres/radiomics/tree/master/STUDIES/LGG_study/WORKSPACE/masterScript_LGG.m>
% for a complete example of how to utilize the current function.
% -------------------------------------------------------------------------
% OUTPUTS: Feature sets are saved in a folder named 'FSET'in the
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

% INITIALIZATION
time = 60; % Number of seconds to wait before checking if parallel computations are done
cd(pathExperiments), load('training'), mkdir('FSET'), cd('FSET'), pathSet = pwd; 
mkdir('batchLog_FSET'), cd('batchLog_FSET'), pathBatch = pwd;
nameOutcomes = fieldnames(training); nOutcomes = numel(nameOutcomes);
for o = 1:nOutcomes
    outcomes.(nameOutcomes{o}) = training.(nameOutcomes{o}).outcome;
end
setNames = fieldnames(training.(nameOutcomes{1}).text);
[param] = batchExperiments(setNames,outcomes,nBatch); nBatch = length(param);

% PRODUCE BATCH COMPUTATIONS
save('workspace','pathSet','pathMINE','training','param','setSize','alpha','delta','nBoot','seed'), pause(5);
for i = 1:nBatch
    nameScript = ['batch',num2str(i),'_script.m'];
    fid = fopen(nameScript,'w');
    fprintf(fid,'load(''workspace'')\n');
    for j = 1:numel(param{i})
        fprintf(fid,['calcAllFeatureSets_LGG(pathSet,training,pathMINE,param{',num2str(i),'}{',num2str(j),'},setSize,alpha,delta,nBoot,seed,',num2str(i),')\n']);
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