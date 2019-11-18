function chooseBestModelsCT_HN(pathWORK,fSetName,nameOutcomes,maxOrder)
% -------------------------------------------------------------------------
% function chooseBestModelsCT_HN(pathWORK,fSetName,nameOutcomes,maxOrder)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function requires user input to choose the set types and model orders
% providing the best parsimonious models for all outcomes analyzed in the
% HN study. See ref.[1] for more details. Goal: comparison with the 
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
% - fSetName: Cell of strings specifying the name of the type of feature 
%             set analyzed (e.g., {'PET', 'SEPARATE', 'FUSED'})
% - nameOutcomes: Cell of strings specifying the outcome names to analyze.
%                 (e.g., {'Failure','Locoregional','Distant','Death'})
% - maxOrder: Integer specifying the maximal multivariable model order.
% -------------------------------------------------------------------------
% OUTPUTS: Final prediction models saved in a folder named 'FINAL_MODELS'
%          in the HN WORKSPACE.
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
cd(pathWORK), mkdir('FINAL_MODELS')
nOutcomes = length(nameOutcomes);
nType = length(fSetName);

for o = 1:nOutcomes
    pathOutcome = [pathWORK,'/RESULTS/',nameOutcomes{o}];
    fprintf('\n*********************** DISPLAYING PREDICTION RESULTS FOR ''%s'' OUTCOME ***********************\n',upper(nameOutcomes{o}))
    plotPredictionResults_HN(pathOutcome,fSetName,{'AUC632','Sensitivity632','Specificity632'},maxOrder)
    if nType > 1
        while 1
            setString = ['Which feature set provides the best parsimonious model? \n'];
            for i = 1:nType
                setString = [setString,'--> For the ''',fSetName{i},''' feature set, type ''',num2str(i),''' and press ENTER \n'];
            end
            setString = [setString,'ANSWER: '];
            set = input(setString);
            fprintf('\n')
            if isnumeric(set)
                set = floor(set);
                if set <= nType && set > 0
                    break
                end
            end
        end
    else
        set = 1;
    end
    while 1
        order = input(['Which model order of the ',fSetName{set},' feature set provides the best parsimonious model? \n' ...
                       '--> Type a number between 1 to ',num2str(maxOrder),' and press ENTER \n' ...
                       'ANSWER: ']);
        fprintf('\n')
        if isnumeric(order) && order <= 10 && order >= 1
            break
        end
    end
    cd([pathWORK,'/RESULTS/',nameOutcomes{o}])
    results = load(['RESULTS_',fSetName{set},'_BEST']); results = struct2cell(results); results = results{1}; 
    finalModel = results.(['Order',num2str(order)]);
    finalModel.outcome = nameOutcomes{o};
    cd([pathWORK,'/FINAL_MODELS/',nameOutcomes{o}]), mkdir('CTonly'), cd('CTonly'), save('finalModel','finalModel')
    close all
end

cd(startpath)
end