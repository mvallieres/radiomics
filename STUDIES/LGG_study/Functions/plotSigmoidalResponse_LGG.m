function plotSigmoidalResponse_LGG(resp,Y,modelCI,nameOutcome,pathFig)
% -------------------------------------------------------------------------
% function plotSigmoidalResponse_LGG(resp,Y,modelCI,nameOutcome,pathFig)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function plots the sigmoidal response (P(y=1|x)) of a linear 
% multivariable model response. See ref. [1] for more details.
% -------------------------------------------------------------------------
% REFERENCE:
% [1] Vallieres, M. et al. (2015). A radiomics model from joint FDG-PET and 
%     MRI texture features for the prediction of lung metastases in soft-tissue 
%     sarcomas of the extremities. Physics in Medicine and Biology, 60(14), 
%     5471-5496. doi:10.1088/0031-9155/60/14/5471
% -------------------------------------------------------------------------
% INPUTS:
% 1. resp: Column vector of size [nInst X 1] specifying the linear 
%          multivariable model response, where 'nInst' refers to the number
%          of instances.
% 2. Y: Column vector of size [nInst X 1] specifying the outcome status 
%       (1 or 0) for all instances.
% 3. modelCI: Column vector of size [nInst X 2] specifying the 95%
%             confidence interval on the multivariable model response as
%             defined by the 2.5 (modelCI(i,1)) and the 97.5 (modelCI(i,2))
%             percentiles, for the ith instance. To plot without CIs, use []
%             for the 3rd argument.See ref. [1] for more 
%             details.
% 4. nameOutcome: String specifying the name of the
%                 modeled outcome.
%                 --> Ex: 'progression'
% 5. pathFig: (optional).  Full path to where figure is saved without
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

startpath = pwd;

if nargin < 5
    pathFig = '';
end

prob = 1./(1 + exp(-resp)); % Sigmoidal response, or P(y=1|x)
respPos = resp(Y==1); 
respNeg = resp(Y==0); 
probPos = prob(Y==1);
probNeg = prob(Y==0);
name = nameOutcome;
symbols = {'ob','xr'}; % First entry for positive instances, second entry for negative instances

if isempty(pathFig)
    figure
else
    h = figure('visible','off');
end
[sortResp,~] = sort(resp);
sigX = (sortResp(1)-5):0.1:sortResp(end)+5;
sigY = 1./(1 + exp(-sigX));
plot(sigX,sigY,'-k','LineWidth',3);
hold on
if ~isempty(modelCI)
    lowXpos = modelCI(Y==1,1); highXpos = modelCI(Y==1,2);
    lowXneg = modelCI(Y==0,1); highXneg = modelCI(Y==0,2);
    [~] = herrorbar(respPos,probPos,respPos-lowXpos,highXpos-respPos,'ob');
    hold on
    [~] = herrorbar(respNeg,probNeg,respNeg-lowXneg,highXneg-respNeg,'xr');
    legend('Sigmoidal response',['95% CI: ',name,' = 1'],['Status: ',name,' = 1'],['95% CI: ',name,' = 0'],['Status: ',name,' = 0'],'Location','NorthWest')
else
    plot(respPos,probPos,symbols{1},'LineWidth',6,'MarkerSize',18,'MarkerFaceColor',symbols{1}(end),'MarkerEdgeColor',symbols{1}(end));
    hold on
    plot(respNeg,probNeg,symbols{2},'LineWidth',6,'MarkerSize',20,'MarkerFaceColor',symbols{2}(end),'MarkerEdgeColor',symbols{2}(end));
    legend('Sigmoidal response',['Status: ',name,' = 1'],['Status: ',name,' = 0'],'Location','NorthWest')
end
xlabel('Multivariable model response','FontSize',30)
ylabel(['Probability that ',name,' = 1'],'FontSize',30)
axis([sortResp(1)-0.33*abs(sortResp(1)),sortResp(end)+0.33*abs(sortResp(end)),-0.05,1.05])
set(gca,'FontSize',24)

if ~isempty(pathFig)
    cd(pathFig)
    saveas(h,['sigPlotResponse_',name],'fig')
end

cd(startpath)
end
