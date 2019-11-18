function plotPredictionResults_HN_tempFig3(pathResults,nameOutcome,fSetNames,metric,maxOrder)
% -------------------------------------------------------------------------
% function plotPredictionResults_HN(pathResults,nameOutcome,fSetName,metrics,maxOrder)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function plots prediction performance estimation results for all the
% different feature set types entered as inputs. See ref. [1] for more
% details.
% -------------------------------------------------------------------------
% REFERENCE:
% [1] Vallieres, M. et al. (2015). FDG-PET/CT radiomics models for the 
%     early prediction of different tumour outcomes in head and neck cancer.
%     The Journal of Nuclear Medicine, aa(bb), xxx-yyy. 
%     doi:
% -------------------------------------------------------------------------
% INPUTS:
% 1. pathResults: Full path to the 'RESULTS' folder where prediction results
%                 are saved.
%                 --> Ex: '/myProject/WORKSPACE/CV-BASED_RESULTS/Experiment1/RESULTS'
% 2. nameOutcome: String specifying the name of the outcome being displayed
%                 --> Ex: 'Distant'
% 3. fSetName: String specifying the feature set name being displayed
%              --> Ex: 'FUSED'
% 4. metrics: Cell of strings specifying the names of the metrics to show 
%             in the plot. 
%             --> Ex: {'AUC632','Sensitivity632','Specificity632'}
% 5. maxOrder: Integer specifying the maximal multivariable model order.
%              --> Ex: 10
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


% THIS FUNCTION FINALLY PLOTS ONLY THE AUC OR CI, ADAPT HEADER FOR THIS

startpath = pwd;
cd(pathResults)

signs = {'-r',':b','-.g','-.y'};
if iscell(fSetNames)
    nFSET = numel(fSetNames);
else
    nFSET = 1;
end

maxOrderChosen = maxOrder;


for i = 1:nFSET
    if nFSET > 1
        fSET = fSetNames{i};
    else
        fSET = fSetNames;
    end
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
    errorbar(1:maxOrder,val(:,1),val_SE(:,1),signs{i},'LineWidth',4,'MarkerFaceColor',signs{i}(end),'MarkerSize',8)
    hold on
end
set(gca,'FontSize',20)
xlabel('Model Order','FontSize',26)
%ylabel('Prediction estimation','FontSize',20)
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
ylabel(metric,'FontSize',26)
%title([nameOutcome,' -- ',metric],'FontSize',24,'FontWeight','bold')
title([nameOutcome],'FontSize',30,'FontWeight','bold')
legend(fSetNames,'Location','SouthEast')
axis([0 maxOrderChosen+1 0.5 1])
set(gca,'XTick',[1 2 3 4 5 6 7 8 9 10])
axis square, grid on

cd(startpath)
end