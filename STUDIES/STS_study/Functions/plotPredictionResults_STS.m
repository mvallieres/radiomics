function plotPredictionResults_STS(pathResults,fSetName,metrics,maxOrder)
% -------------------------------------------------------------------------
% function plotPredictionResults_STS(pathResults,fSetName,metrics,maxOrder)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function plots prediction performance estimation results for all the
% different feature set types entered as inputs. See ref. [1] for more
% details.
% -------------------------------------------------------------------------
% REFERENCE:
% [1] Vallieres, M. et al. (2015). A radiomics model from joint FDG-PET and 
%     MRI texture features for the prediction of lung metastases in soft-tissue 
%     sarcomas of the extremities. Physics in Medicine and Biology, 60(14), 
%     5471-5496. doi:10.1088/0031-9155/60/14/5471
% -------------------------------------------------------------------------
% INPUTS:
% - pathResults: Full path to the 'RESULTS' folder where prediction results
%                are saved.
% - fSetName: Cell of strings specifying the name of the type of feature 
%                 sets to analyze.
%                 Example: fSetName = {'PET','SEPARATE','FUSED'}
% - metrics: Name of the metrics to show in the plots. 
%            Example: metrics = {'AUC632','Sensitivity632','Specificity632'}
% - maxOrder: Integer specifying the maximal multivariable model order. 
% -------------------------------------------------------------------------
% AUTHOR(S): Martin Vallieres <mart.vallieres@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: May 2015
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
cd(pathResults)

signs = {'-r',':b','--g'};
nFSET = numel(fSetName);
nMetrics = numel(metrics);
figure
for i = 1:nFSET
    val = zeros(maxOrder,numel(metrics));
    val_SE = zeros(maxOrder,nMetrics);
    results = load(['RESULTS_',fSetName{i},'_BEST']); results = struct2cell(results); results = results{1};
    for j = 1:maxOrder
        orderName = ['Order',num2str(j)];
        for k =1:nMetrics
            val(j,k) = results.(orderName).(metrics{k});
            val_SE(j,k) = results.(orderName).(['SE_',metrics{k}]);
        end
    end
    subplot(1,nFSET,i)
    for k = 1:nMetrics
        errorbar(1:maxOrder,val(:,k),val_SE(:,k),signs{k},'LineWidth',2,'MarkerFaceColor',signs{k}(end),'MarkerSize',4)
        hold on
    end
    xlabel('Model Order','FontSize',24)
    ylabel('Prediction performance','FontSize',24)
    title(fSetName{i},'FontSize',30,'FontWeight','bold')
    legend(metrics,'Location','SouthEast')
    axis([0 maxOrder+1 0.5 1])
    set(gca,'FontSize',20)
    set(gca,'XTick',[1 2 3 4 5 6 7 8 9 10])
    hold off
end

cd(startpath)
end