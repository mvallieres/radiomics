function plotSigmoidalResponse(resp,Y,modelCI,nameOutcome)
% -------------------------------------------------------------------------
% function plotSigmoidalResponse(resp,Y,modelCI,nameOutcome)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function plots the sigmoidal response (P(y=1|x)) of a linear 
% multivariable model response. See ref. [1] for more details.
% -------------------------------------------------------------------------
% REFERENCE:
% [1] Vallieres, M. et al. (2015). A radiomics model from joint FDG-PET and 
%     MRI texture features for the prediction of lung metastases in soft-tissue 
%     sarcomas of the extremities. Physics in Medicine and Biology, 00(0), xxx-yyy. 
% -------------------------------------------------------------------------
% INPUTS:
% - resp: Column vector of size [nInst X 1] specifying the linear 
%         multivariable model response, where 'nInst' refers to the number
%         of instances.
% - Y: Column vector of size [nInst X 1] specifying the outcome status 
%      (1 or 0) for all instances.
% - modelCI: Column vector of size [nInst X 2] specifying the 95%
%            confidence interval on the multivariable model response as
%            defined by the 2.5 (modelCI(i,1)) and the 97.5 (modelCI(i,2))
%            percentiles, for the ith instance. See ref. [1] for more 
%            details.
% - nameOutcome: (optional input). String specifying the name of the
%                modeled outcome.
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

prob = 1./(1 + exp(-resp)); % Sigmoidal response, or P(y=1|x)
respPos = resp(Y==1); 
respNeg = resp(Y==0); 
probPos = prob(Y==1);
probNeg = prob(Y==0);
lowXpos = modelCI(Y==1,1); highXpos = modelCI(Y==1,2);
lowXneg = modelCI(Y==0,1); highXneg = modelCI(Y==0,2);
if nargin == 4
    name = nameOutcome;
else
    name = 'Y';
end


figure
[sortResp,ind] = sort(resp);
h = plot(sortResp,prob(ind),'-k','LineWidth',3);
hold on
h = herrorbar(respPos,probPos,lowXpos,highXpos,'ob');
hold on
h = herrorbar(respNeg,probNeg,lowXneg,highXneg,'xr');
xlabel('Multivariable model response','FontSize',24)
ylabel(['Probability that ',name,' = 1'],'FontSize',24)
legend('Sigmoidal response',['95% CI: ',name,' = 1'],['Status: ',name,' = 1'],['95% CI: ',name,' = 0'],['Status: ',name,' = 0'],'Location','NorthWest')
axis([sortResp(1)-5,sortResp(end)+5,-0.05,1.05])
set(gca,'FontSize',20)

end
