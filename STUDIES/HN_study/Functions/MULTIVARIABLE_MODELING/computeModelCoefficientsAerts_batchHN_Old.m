function computeModelCoefficientsAerts_batchHN_Old(pathExperiments,nExp,imbalance,nBatch,matlabPATH,seed)
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
time = 30; % Number of seconds to wait before checking if parallel computations are done

for exp = 1:nExp
    cd(fullfile(pathExperiments,['Experiment',num2str(exp)])), load('training')
    mkdir('AertsSign'), cd('AertsSign'), pathFinalModels = pwd;
    mkdir('batchLog_Coeff'), cd('batchLog_Coeff'), pathBatch = pwd; cd(pathFinalModels)
    nameOutcomes = fieldnames(training); nOutcomes = numel(nameOutcomes);
    setNames = fieldnames(training.(nameOutcomes{1}).sign); nFset = numel(setNames);
    for o = 1:nOutcomes
        outcomes.(nameOutcomes{o}) = training.(nameOutcomes{o}).outcome;
        mkdir(nameOutcomes{o}), cd(nameOutcomes{o})
        for f = 1:nFset
            mkdir(setNames{f}), cd(setNames{f})
            finalModel = struct;
            finalModel.Data = training.(nameOutcomes{o}).sign.(setNames{f});
            finalModel.Name = {'Energy';'Compactness';'GLN';'GLN_HLH'};
            finalModel.Order = 4;
            finalModel.outcome = nameOutcomes{o};
            save('finalModel','finalModel')
            cd ..
        end
        cd ..
    end
    [param] = batchExperiments(setNames,outcomes,nBatch); nBatch = length(param);
    cd(pathBatch)
    
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
            fprintf(fid,['[coeff,response,modelCI] = computeModelCoefficients_HN(finalModel.Data,training.',param{i}{j}{2},'.outcome,imbalance,seed);\n']);
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
end

cd(startpath)
end