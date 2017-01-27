function prepareFsetVASARI(pathWORK,pathSTUDY)
% -------------------------------------------------------------------------
% function prepareFsetVASARI(pathWORK,pathSTUDY)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function prepares feature sets for multivariable modeling
% experiments with VASARI features
% -------------------------------------------------------------------------
% INPUTS:
% 1. pathWORK: Full path to the WORKSPACE of the LGG study. Defined in
%              masterScript_LGG.m. 
%              --> Ex: '/myProject/WORKSPACE/'
% 2. pathSTUDY: Full path to where the organized data will be saved.
%               Defined in masterScript_LGG.m
%               --> Ex: '/myProject/WORKSPACE/STUDY_DATA'
% -------------------------------------------------------------------------
% OUTPUTS: Organized 'feature sets data in
% '/myProject/WORKSPACE/VASARI/FSET
% -------------------------------------------------------------------------
% AUTHOR(S): 
% - Martin Vallieres <mart.vallieres@gmail.com>
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

cd(pathSTUDY), load('outcomes')
nameOutcomes = fieldnames(outcomes); nOutcomes = numel(nameOutcomes);

cd(pathWORK), mkdir('VASARI'), cd('VASARI'), save('outcomes','outcomes'), mkdir('FSET'), cd('FSET'), pathFset = pwd;

for o = 1:nOutcomes
    outcomeName = nameOutcomes{o};
    cd(pathSTUDY), load(['vasari_',outcomeName]) % Variable 'vasari' now in MATLAB workspace
    outcome = outcomes.(outcomeName);
    vasari(isnan(vasari)) = 0; nVas = size(vasari,2); nameVas = cell(nVas,1);
    for i = 1:nVas
        nameVas{i} = ['F',num2str(i)];
    end
    delete = [];
    for i = 1:nVas
        if numel(unique(vasari(:,i))) < 3 % Deleted, not enough variation
            delete = [delete,i];
        end
    end
    vasari(:,delete) = []; nameVas(delete) = [];
    
    fSet = struct;
    fSet.Data = vasari;
    fSet.Info = nameVas;
    cd(pathFset), save(['FSET_VASARI_',outcomeName],'fSet')
end

cd(startpath)
end