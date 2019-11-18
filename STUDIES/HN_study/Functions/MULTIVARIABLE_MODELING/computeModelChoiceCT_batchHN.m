function computeModelChoiceCT_batchHN(pathWORK,fSetName,outcomes,freedomMat,maxOrder,nBoot,imbalance,nBatch,matlabPATH)
% -------------------------------------------------------------------------
% function computeModelChoiceCT_batchHN(pathWORK,fSetName,outcomes,freedomMat,maxOrder,nBoot,imbalance,nBatch,matlabPATH)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function computes feature set selection for the 'CT' feature set
% only and for all outcomes. See ref. [1,2] for more details. Goal: comparison 
% with the "Radiomics signature" of (Aerts et al., Nat Commun, 2014).
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
% - fSetName: Cell of strings specifying the name of the type of feature 
%             set analyzed (e.g., {'PET', 'SEPARATE', 'FUSED'})
% - outcomes: Structure specifying the status (1 or 0) for different
%             outcomes in HN cancer. Contains: outcomes.Failure, 
%             outcomes.Locoregional, outcomes.Distant, outcomes.Death. See
%             ref.[1] for more details.
% - freedomMat:  Cell vector of matrices of row vectors of 1's and 0's to 
%                specify the degree of freedom on texture extraction 
%                parameters for all experiments. For example, for an nth 
%                experiment of an ith scan where extraction parameters 1, 2
%                and 4 in paramCell are allowed to vary, use 
%                freedomMat{i}(n,:) = [1,1,0,1].
% - maxOrder: Integer specifying the maximal model order to construct.
% - nBoot: Number of bootstrap samples to use.
% - imbalance: String specifying the type of imbalance-adjustement strategy
%              employed. Either 'IABR' for imbalance-adjusted bootstrap
%              resampling (see ref.[1]), or 'IALR' for imbalance-adjusted
%              logistic regression (see ref.[2]).
% - nBatch: Number of parallel batch.
% - matlabPATH: Full path to the MATLAB excutable on the system.
% -------------------------------------------------------------------------
% OUTPUTS: Mutivariable models are saved in a folder named 'MODELS' in the 
%          HN WORKSPACE.
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

% INITIALIZATON
cd(pathWORK), cd('MODELS')
nameOutcomes = fieldnames(outcomes);
nOutcomes = length(nameOutcomes);
outcomeSep = cell(1,nOutcomes);
for o = 1:nOutcomes
    temp.(nameOutcomes{o}) = outcomes.(nameOutcomes{o});
    outcomeSep{o} = temp;
    temp = [];
end
nType = length(fSetName);
[param] = batchExperiments(nType,outcomes,nBatch); nBatch = length(param);
cd('SCRIPTS'), pathScripts = pwd;
time = 60; % Number of seconds to wait before checking if parallel computations are done


% PRODUCE BATCH COMPUTATIONS
save('workspace'), pause(5);
for i = 1:nBatch
    nameScript = ['Script_FeatSel_CT',num2str(i),'.m'];
    fid = fopen(nameScript,'w');
    fprintf(fid,'load(''workspace'')\n');
    for j = 1:size(param{i},1)
        fprintf(fid,['computeAllModelChoice_HN(pathWORK,fSetName{param{',num2str(i),'}(',num2str(j),',1)},outcomeSep{param{',num2str(i),'}(',num2str(j),',2)},freedomMat{param{',num2str(i),'}(',num2str(j),',1)},maxOrder,nBoot,imbalance,',num2str(i),')\n']);
    end
    fprintf(fid,'clear all');
    fclose(fid);
    system([matlabPATH,' -nojvm -nodisplay -nodesktop -nosplash < ',nameScript,' >& ',nameScript(1:end-1),'log &']);
end


% WAITING LOOP
nameVerif = {};
for i = 1:nType
    nExp = size(freedomMat{i},1);
    nParam = size(freedomMat{i},2);
    for j = 1:nExp
        temp = ['MODELS_',fSetName{i},'_'];
        for k = 1:nParam
            temp = [temp,num2str(freedomMat{i}(j,k))];
        end
        nameVerif = [nameVerif;temp];
    end
end
nVerif = numel(nameVerif);
while 1
    pause(time);
    check = zeros(nVerif,1);
    for o = 1:nOutcomes
        cd([pathWORK,'/MODELS/',nameOutcomes{o}])
        for i = 1:nVerif
            check(i) = check(i) + exist([nameVerif{i},'.mat']);
        end
    end
    if sum(check) == nVerif*nOutcomes*2
        break
    end
end

cd(pathScripts)
delete('workspace.mat')
cd(startpath)
end