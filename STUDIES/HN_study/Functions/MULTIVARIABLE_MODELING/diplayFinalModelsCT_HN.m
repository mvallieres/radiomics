function diplayFinalModelsCT_HN(pathWORK,outcomes)
% -------------------------------------------------------------------------
% function diplayFinalModelsCT_HN(pathWORK,outcomes)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function displays the final models for all analyzed outcomes in the
% HN study for the 'CT' feature set only, in terms of: multivariable model 
% response, model variables, and prediction performance estimation. See 
% ref.[1] for more details. Goal: comparison with the "Radiomics signature" 
% of (Aerts et al., Nat Commun, 2014)
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
nameOutcomes = fieldnames(outcomes);
nOutcomes = length(nameOutcomes);

for o = 1:nOutcomes
    cd([pathWORK,'/FINAL_MODELS/',nameOutcomes{o},'/CTonly'])
    load('finalModel'), load('coeff'), load('modelCI'), load('response') % Variables 'finalModel', 'coeff',  'modelCI', and 'response' now present in MATLAB workspace
    order = size(finalModel.Data,2);
    plotSigmoidalResponse(response,outcomes.(nameOutcomes{o}),modelCI,upper(nameOutcomes{o}))
    fprintf(['\n\n --> THE FINAL MULTIVARIABLE MODEL (CT FEATURE SET) FOR ''%s'' OUTCOME IS:\n\n'...
             '               g(x) =               \n'],upper(nameOutcomes{o}))
    for i = 1:order
        fprintf([num2str(coeff(i)),' X ',finalModel.Name{i},'\n'])
        fprintf('                    +               \n')
    end
    fprintf(['                   ',num2str(coeff(end)),'\n'])
    fprintf('\nWITH CORRESPONDING PREDICTION PERFORMANCE ESTIMATION:\n')
    fprintf(['AUC = ',num2str(roundsd(finalModel.AUC632,ceil(log10(finalModel.AUC632/roundsd(finalModel.SE_AUC632,1))))),' ± ',num2str(roundsd(finalModel.SE_AUC632,1)),'\n'])
    fprintf(['Sensitivity = ',num2str(roundsd(finalModel.Sensitivity632,ceil(log10(finalModel.Sensitivity632/roundsd(finalModel.SE_Sensitivity632,1))))),' ± ',num2str(roundsd(finalModel.SE_Sensitivity632,1)),'\n'])
    fprintf(['Specificity = ',num2str(roundsd(finalModel.Specificity632,ceil(log10(finalModel.Specificity632/roundsd(finalModel.SE_Specificity632,1))))),' ± ',num2str(roundsd(finalModel.SE_Specificity632,1)),'\n'])
    fprintf('\n\n')
end

cd(startpath)
end