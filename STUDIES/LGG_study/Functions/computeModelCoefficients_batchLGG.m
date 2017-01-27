function computeModelCoefficients_batchLGG(pathExperiments,imbalance,nBatch,matlabPATH,seed)
% -------------------------------------------------------------------------
% function computeModelCoefficients_batchLGG(pathExperiments,imbalance,nBatch,matlabPATH,seed)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function computes the final logistic regression coefficients and
% bootstrap confidence intervals of the final models obtained for all
% outcomes analyzed in the LGG study.
% -------------------------------------------------------------------------
% INPUTS:
% 1. pathExperiments: Full path to the directory containing all experiments.
%                     --> Ex: '/myProject/WORKSPACE/LOGISTIC_REGRESSION'
% 2. imbalance: String specifying the type of imbalance-adjustement strategy
%               employed. Either 'IABR' for imbalance-adjusted bootstrap
%               resampling (see ref.[1]), or 'IALR' for imbalance-adjusted
%               logistic regression (formal reference to come).
%               --> Ex: 'IALR'
% 3. nBatch: Number of parallel batch.
%            --> Ex: 8
% 4. matlabPATH: Full path to the MATLAB executable on the system.
%                --> 'matlab' if a symbolic link to the matlab executable
%                     was previously created.
% 5. seed: Numerical number to use as seed for bootstrapping experiment
%          --> Ex: 54288
% -------------------------------------------------------------------------
% OUTPUTS: Final coefficients, model response and confidence intervals 
%          saved in a folder named '/myProject/WORKSPACE/LOGISTIC_REGRESSION/FINAL_MODELS/nameOutcome/fSetName'.
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
time = 30; % Number of seconds to wait before checking if parallel computations are done


cd(pathExperiments), load('training')
cd('FINAL_MODELS'), pathFinalModels = pwd;
mkdir('batchLog_Coeff'), cd('batchLog_Coeff'), pathBatch = pwd;
nameOutcomes = fieldnames(training); nOutcomes = numel(nameOutcomes);
for o = 1:nOutcomes
    outcomes.(nameOutcomes{o}) = training.(nameOutcomes{o}).outcome;
end
setNames = fieldnames(training.(nameOutcomes{1}).text);
[param] = batchExperiments(setNames,outcomes,nBatch); nBatch = length(param);

% PRODUCE BATCH COMPUTATIONS
save('workspace','pathFinalModels','training','param','pathBatch','imbalance','seed'), pause(3)
for i = 1:nBatch
    nameScript = ['batch',num2str(i),'_script.m'];
    fid = fopen(nameScript,'w');
    fprintf(fid,'tic\n');
    fprintf(fid,'load(''workspace'')\n');
    for j = 1:numel(param{i})
        fprintf(fid,['cd(fullfile(pathFinalModels,param{',num2str(i),'}{',num2str(j),'}{2},param{',num2str(i),'}{',num2str(j),'}{1}))\n']);
        fprintf(fid,'load(''finalModel'')\n');
        fprintf(fid,['fprintf(''COMPUTING THE LOGISTIC REGRESSION COEFFICIENTS OF THE FINAL MODEL OF "',param{i}{j}{2},'" OUTCOME, "',param{i}{j}{1},'" FEATURE SET ... '')']);
        fprintf(fid,'\n');
        fprintf(fid,['[coeff,response,modelCI] = computeModelCoefficients(finalModel.Data,training.',param{i}{j}{2},'.outcome,imbalance,seed);\n']);
        fprintf(fid,['fprintf(''DONE!\\n'')']);
        fprintf(fid,'\n');
        fprintf(fid,'save(''coeff'',''coeff''), save(''response'',''response''), save(''modelCI'',''modelCI'')\n');
    end
    fprintf(fid,'cd(pathBatch)\n');
    fprintf(fid,['system(''touch batch',num2str(i),'_end'');\n']);
    fprintf(fid,'clear all\n');
    fprintf(fid,'toc');
    fclose(fid);
    system([matlabPATH,' -nojvm -nodisplay -nodesktop -nosplash < ',nameScript,' >& ',nameScript(1:end-1),'log &']);
end
waitBatch(pathBatch,time,nBatch)
delete('workspace.mat')


cd(startpath)
end