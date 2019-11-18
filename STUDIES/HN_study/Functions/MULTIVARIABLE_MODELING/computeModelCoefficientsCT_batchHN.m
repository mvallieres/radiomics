function computeModelCoefficientsCT_batchHN(pathWORK,outcomes,imbalance,matlabPATH)
% -------------------------------------------------------------------------
% function computeModelCoefficientsCT_batchHN(pathWORK,outcomes,imbalance,matlabPATH)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function computes the final logistic regression coefficients and
% bootstrap confidence intervals of the final models obtained for all
% outcomes analyzed in the HN study for the 'CT' feature set only. 
% See ref.[1] for more details. Goal: comparison with the 
% "Radiomics signature" of (Aerts et al., Nat Commun, 2014)
% -------------------------------------------------------------------------
% REFERENCE:
% [1] Vallieres, M. et al. (2015). FDG-PET/CT radiomics models for the 
%     early prediction of different tumour outcomes in head and neck cancer.
%     The Journal of Nuclear Medicine, aa(bb), xxx-yyy. 
%     doi:
% -------------------------------------------------------------------------
% INPUTS:
% - pathWORK: Full path to the HN WORKSPACE directory.
% - outcomes: Structure specifying the status (1 or 0) for different
%             outcomes in HN cancer. Contains: outcomes.Failure, 
%             outcomes.Locoregional, outcomes.Distant, outcomes.Death. See
%             ref.[1] for more details.
% - imbalance: String specifying the type of imbalance-adjustement strategy
%              employed. Either 'IABR' for imbalance-adjusted bootstrap
%              resampling, or 'IALR' for imbalance-adjusted logistic
%              regression.
% - matlabPATH: Full path to the MATLAB excutable on the system.
% -------------------------------------------------------------------------
% OUTPUTS: Final coefficients, model response and confidence intervals 
%          saved in a folder named 'FINAL_MODELS' in the HN WORKSPACE.
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
cd(pathWORK)
nameOutcomes = fieldnames(outcomes);
nOutcomes = length(nameOutcomes);
time = 60; % Number of seconds to wait before checking if parallel computations are done

% PRODUCE BATCH COMPUTATIONS
for o = 1:nOutcomes
    cd([pathWORK,'/FINAL_MODELS/',nameOutcomes{o},'/CTonly'])
    save('workspace'), pause(5);
    nameScript = ['Script_CoeffCT_',nameOutcomes{o},'.m'];
    fid = fopen(nameScript,'w');
    fprintf(fid,'tic\n');
    fprintf(fid,'load(''workspace'')\n');
    fprintf(fid,'load(''finalModel'')\n');
    fprintf(fid,['fprintf(''COMPUTING THE LOGISTIC REGRESSION COEFFICIENTS OF THE FINAL MODEL OF ''''',upper(nameOutcomes{o}),''''' OUTCOME ... '')']);
    fprintf(fid,'\n');
    fprintf(fid,['[coeff,response,modelCI] = computeModelCoefficients(finalModel.Data,outcomes.(nameOutcomes{',num2str(o),'}),imbalance,',num2str(o),');\n']);
    fprintf(fid,['fprintf(''DONE!\\n'')']);
    fprintf(fid,'\n');
    fprintf(fid,'save(''coeff'',''coeff''), save(''response'',''response''), save(''modelCI'',''modelCI'')\n');
    fprintf(fid,'clear all\n');
    fprintf(fid,'toc');
    fclose(fid);
    system([matlabPATH,' -nojvm -nodisplay -nodesktop -nosplash < ',nameScript,' >& ',nameScript(1:end-1),'log &']);
end

% WAITING LOOP
while 1
    pause(time);
    check = zeros(nOutcomes,1);
    for o = 1:nOutcomes
        cd([pathWORK,'/FINAL_MODELS/',nameOutcomes{o},'/CTonly'])
        check(o) = exist('modelCI.mat');
    end
    if sum(check) == nOutcomes*2
        break
    end
end

for o = 1:nOutcomes
    cd([pathWORK,'/FINAL_MODELS/',nameOutcomes{o},'/CTonly'])
    delete('workspace.mat')
end

cd(startpath)
end