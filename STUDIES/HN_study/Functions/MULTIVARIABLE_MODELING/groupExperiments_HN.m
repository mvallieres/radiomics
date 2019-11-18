function groupExperiments_HN(pathWORK,nameOutcomes,fSetName,freedomMat,maxOrder)
% -------------------------------------------------------------------------
% function groupExperiments_HN(pathWORK,nameOutcomes,fSetName,freedomMat,maxOrder)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function finds the best combination of model order and texture 
% extraction parameter degree of freedom in terms of [AUC]0.632+ for all 
% experiments performed for all feature sets. See ref. [1,2] for more 
% details.
% -------------------------------------------------------------------------
% REFERENCES:
% [1] Vallieres, M. et al. (2015). FDG-PET/CT radiomics models for the 
%     early prediction of different tumour outcomes in head and neck cancer.
%     The Journal of Nuclear Medicine, aa(bb), xxx-yyy. 
%     doi:
% [2] Vallieres, M. et al. (2015). A radiomics model from joint FDG-PET and 
%     MRI texture features for the prediction of lung metastases in soft-tissue 
%     sarcomas of the extremities. Physics in Medicine and Biology, 60(14), 
%     5471-5496. doi:10.1088/0031-9155/60/14/5471
% -------------------------------------------------------------------------
% INPUTS:
% - pathWORK: Full path to the HN WORKSPACE directory.
% - nameOutcomes: Cell of strings specifying the the names of the outcomes
%                 modeled in the HN study. 
%                 Example: {'Failure','Locoregional','Distant','Death'}
% - fSetName: Cell of strings specifying the name of the type of feature 
%             set analyzed (e.g., 'PET', 'SEPARATE', 'FUSED', etc.)
% - freedomMat:  Cell vector of matrices of row vectors of 1's and 0's to 
%                specify the degree of freedom on texture extraction 
%                parameters for all experiments. For example, for an nth 
%                experiment of an ith scan where extraction parameters 1, 2
%                and 4 in paramCell are allowed to vary, use 
%                freedomMat{i}(n,:) = [1,1,0,1].
% - maxOrder: Integer specifying the maximal model order to construct.
% -------------------------------------------------------------------------
% OUTPUTS: The optimal prediction performance results for every model order
%          (in terms of etexture extraction parameter degree of freedom) 
%          are saved in a folder named 'RESULTS' in the HN WORKSPACE.
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
nOutcomes = length(nameOutcomes);
nType = length(fSetName);

for o = 1:nOutcomes
    cd([pathWORK,'/RESULTS/',nameOutcomes{o}])
    for n = 1:nType
        
        % LOADING ALL RESULTS FOR THE DIFFERENT EXPERIMENTS
        nFreedom = size(freedomMat{n},1);
        nParamType = size(freedomMat{n},2);
        cellResults = cell(1,nFreedom);
        for i = 1:nFreedom
            nameOpen = ['RESULTS_',fSetName{n},'_'];
            for j = 1:nParamType
                nameOpen = [nameOpen,num2str(freedomMat{n}(i,j))];
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
            results.(orderName).freedom = freedomMat{n}(indMax,:);
        end

        % SAVING RESULTS
        save(['RESULTS_',fSetName{n},'_BEST'],'results')
    
    end
end

cd(startpath)
end