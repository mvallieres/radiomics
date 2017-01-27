function plotPredictionResults_LGG(pathResults,nameOutcome,fSetNames,metric,maxOrder,pathFig)
% -------------------------------------------------------------------------
% function plotPredictionResults_LGG(pathResults,nameOutcome,fSetNames,metrics,maxOrder)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function plots prediction performance estimation results for all the
% different feature set types entered as inputs in the LGG study.
% -------------------------------------------------------------------------
% INPUTS:
% 1. pathResults: Full path to the 'RESULTS' folder where prediction results
%                 are saved.
%                 --> Ex: '/myProject/WORKSPACE/LOGISTIC_REGRESSION/RESULTS'
% 2. nameOutcome: String specifying the name of the outcome being displayed
%                 --> Ex: 'progression'
% 3. fSetNames: Cell of strings specifying the name of the type of feature 
%               set analyzed.
%               --> Ex: {'T1W_T2W','T1W_T2F','T1CE_T2W','T1CE_T2F'}
% 4. metric: String specifying the metric to display.
%            --> 'AUC632'
% 5. maxOrder: Integer specifying the maximal multivariable model order.
%              --> Ex: 10
% 6. pathFig: (optional).  Full path to where figure is saved without
%             displaying it. Put '' for displaying the figure and not 
%             saving it to 'pathFig' (default).
%             --> Ex: ''
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

if nargin < 6
    pathFig = '';
end

startpath = pwd;
cd(pathResults)

signs = {'-r',':b','--g','-.y'};
nFSET = numel(fSetNames);

maxOrderChosen = maxOrder;

if isempty(pathFig)
    figure
else
    h = figure('visible','off');
end
for i = 1:nFSET
    fSET = fSetNames{i};
    results = load(['RESULTS_',fSET,'_',nameOutcome]); results = struct2cell(results); results = results{1};
    nOrders = numel(fieldnames(results));
    if nOrders < maxOrderChosen
        maxOrder = nOrders;
    else
        maxOrder = maxOrderChosen;
    end
    val = zeros(maxOrder,1);
    val_SE = zeros(maxOrder,1);
    for j = 1:maxOrder
        orderName = ['Order',num2str(j)];
        val(j,1) = results.(orderName).(metric);
        val_SE(j,1) = results.(orderName).(['SE_',metric]);
    end
    errorbar(1:maxOrder,val(:,1),val_SE(:,1),signs{i},'LineWidth',3,'MarkerFaceColor',signs{i}(end),'MarkerSize',6)
    hold on
end
set(gca,'FontSize',20)
xlabel('Model Order','FontSize',24)
ylabel('Prediction performance','FontSize',24)
ind = strfind(metric,'632');
if ~isempty(ind)
    metric = [metric,'+'];
end
metric = [metric(1:ind-1),'_{',metric(ind:end),'}'];
ind = strfind(nameOutcome,'Death');
if ~isempty(ind)
    nameOutcome(ind:ind+4) = [];
    nameOutcome = ['Survival',nameOutcome];
end
titleName = [nameOutcome,' -- ',metric];
title(titleName,'FontSize',30,'FontWeight','bold')
legend(fSetNames,'Location','SouthEast')
axis([0 maxOrderChosen+1 0.5 1])
set(gca,'XTick',[1 2 3 4 5 6 7 8 9 10])

if ~isempty(pathFig)
    cd(pathFig)
    saveas(h,titleName,'fig')
end

cd(startpath)
end