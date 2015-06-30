function groupExperiments_STS(pathWORK,fSetNameType,freedomMat,maxOrder)
% -------------------------------------------------------------------------
% function groupExperiments_STS(pathWORK,fSetNameType,freedomMat,maxOrder)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function finds the best combination of model order and texture 
% extraction parameter degree of freedom in terms of [AUC]0.632+ for all 
% experiments performed for a given feature set. See ref. [1] for more 
% details.
% -------------------------------------------------------------------------
% REFERENCE:
% [1] Vallieres, M. et al. (2015). A radiomics model from joint FDG-PET and 
%     MRI texture features for the prediction of lung metastases in soft-tissue 
%     sarcomas of the extremities. Physics in Medicine and Biology, 60(14), 
%     5471-5496. doi:10.1088/0031-9155/60/14/5471
% -------------------------------------------------------------------------
% INPUTS:
% - pathWORK: Full path to the STS WORKSPACE directory.
% - fSetNameType: String specifying the name of the type of feature set 
%                 (e.g., 'PET', 'SEPARATE', 'FUSED', etc.)
% - freedomMat:  Matrix of row vectors of 1's and 0's to specify the degree 
%                of freedom on texture extraction parameters for all 
%                experiments. For example, for an ith experiment where 
%                extraction parameters 1, 2 and 4 in paramAll are allowed 
%                to vary, use freedomMat(i,:) = [1,1,0,1].
% - maxOrder: Integer specifying the maximal model order to construct.
% -------------------------------------------------------------------------
% OUTPUTS: The optimal prediction performance results for every model order
%          (in terms of etexture extraction parameter degree of freedom) 
%          are saved in a folder named 'RESULTS' in the STS WORKSPACE.
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
cd([pathWORK,'/RESULTS'])

% LOADING ALL RESULTS FOR THE DIFFERENT EXPERIMENTS
nFreedom = size(freedomMat,1);
nParamType = size(freedomMat,2);
cellResults = cell(1,nFreedom);
for i = 1:nFreedom
    nameOpen = ['RESULTS_',fSetNameType,'_'];
    for j = 1:nParamType
        nameOpen = [nameOpen,num2str(freedomMat(i,j))];
    end
    temp = load(nameOpen); temp = struct2cell(temp); temp = temp{1};
    cellResults{i} = temp;
end

% FINDING THE BEST COMBINATIONS (model order/degree of freedom)
results = struct;
for i = 1:maxOrder
    orderName = ['Order',num2str(i)];
    auc = zeros(nFreedom,1);
    for j = 1:nFreedom
        auc(j) = cellResults{j}.(orderName).AUC632; 
    end
    [~,indMax] = max(auc);
    results.(orderName) = cellResults{indMax}.(orderName);
    results.(orderName).freedom = freedomMat(indMax,:);
end

% SAVING RESULTS
save(['RESULTS_',fSetNameType,'_BEST'],'results')

cd(startpath)
end