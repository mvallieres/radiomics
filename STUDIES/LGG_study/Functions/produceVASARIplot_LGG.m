function produceVASARIplot_LGG(pathResults,nameOutcomes,maxOrder,pathFig)
% -------------------------------------------------------------------------
% function produceVASARIplot_LGG(pathResults,nameOutcomes,maxOrder,pathFig)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function produce the Figure 4 of the LGG study.
% -------------------------------------------------------------------------
% INPUTS:
% 1. pathResults: Full path to the directory containing all experiments.
%                     --> Ex: '/myProject/WORKSPACE/VASARI/RESULTS'
% 2. nameOutcomes: Cell of strings specifying the outcome names to analyze.
%                  --> Ex: {'nonIDH1','IDHcodel','progression','lowGrade'}
% 3. maxOrder: Integer specifying the maximal multivariable model order.
%              --> Ex: 10
% 4. pathFig: (optional).  Full path to where figure is saved without
%             displaying it. Put '' for displaying the figure and not 
%             saving it to 'pathFig' (default).
%             --> Ex: ''
% -------------------------------------------------------------------------
% OUTPUTS: Final prediction models saved in a folder named 
%          '/myProject/WORKSPACE/VASARI/FINAL_MODELS/nameOutcome/fSetName'
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
cd(pathResults)

if nargin < 4
    pathFig = '';
end

signs = {'-r',':b','--g','-.y'};
nOutcomes = numel(nameOutcomes);

if isempty(pathFig)
    figure
else
    h = figure('visible','off');
end
for o = 1:nOutcomes
    nameOutcome = nameOutcomes{o};
    plotAllPredictionResults_VASARI_LGG(pathResults,nameOutcome,{'VASARI'},'AUC632',maxOrder,signs{o})
    hold on
end
legend(nameOutcomes)
title(['VASARI'],'FontSize',30,'FontWeight','bold')

if ~isempty(pathFig)
    cd(pathFig)
    saveas(h,['VASARI_PLOT'],'fig')
end

cd(startpath)
end