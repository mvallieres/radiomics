function calcAllNonTextureFeatures_batchHN_TestWithoutPercentInactive(pathData,pathNonText,namePT,nameCT,nameROI,outcomes,featType,nBatch,matlabPATH)
% -------------------------------------------------------------------------
% function calcAllNonTextureFeatures_batchHN(pathData,pathNonText,namePT,nameCT,nameROI,nBatch,matlabPATH)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% Computes 11 non-texture features for all head and neck (HN) DICOM 
% imaging data downloaded from The Cancer Imaging Archive (TCIA) website 
% at: <http://dx.doi.org/xxxxxxxxxxxxxxxxxxxxxxxxxxxxx>, and organized in a 
% 'DATA' directory using the function readAllDICOM_HN.m. In addition, the 
% Spearman's rank correlation between each non-texture features and 
% different outcomes in HN cancer. Results are saved as a structure named 
% 'nonTextures.mat' in the HN WORKSPACE.
% -------------------------------------------------------------------------
% REFERENCE:
% [1] Vallieres, M. et al. (2015). FDG-PET/CT radiomics models for the 
%     early prediction of different tumour outcomes in head and neck cancer.
%     The Journal of Nuclear Medicine, aa(bb), xxx-yyy. 
%     doi:
% -------------------------------------------------------------------------
% INPUTS:
% 1. pathData: Full path to the HN sData files directory.
%              --> Ex: '/myProject/WORKSPACE/DATA'
% 2. pathNonText: Full path to the HN non texture features directory.
%              --> Ex: '/myProject/WORKSPACE/FEATURES/NON_TEXTURES'
% 3. namePT: Cell of strings of all PET sData files to read
%            --> Ex: {'HGJ_001_PT.PTscan.mat';'HGJ_022_PT.PTscan.mat'}
% 4. namePT: Cell of strings of all CT sData files to read
%            --> Ex: {'HGJ_001_CT.CTscan.mat';'HGJ_022_CT.CTscan.mat'}
% 5. nameROI: Cell of strings specifying the ROI names to analyze for the
%             patients defined by "namePT" and "nameCT"
%             --> Ex: {'GTV';'GTV-P'}
% 6. outcomes: Structure specifying the status (1 or 0) for different
%              outcomes in HN cancer. Contains: outcomes.Failure, 
%              outcomes.Locoregional, outcomes.Distant. See ref.[1] for 
%              more details.
% 7. featType: Either 'GTVp' for primary GTV, or 'GTVtot' for primaty GTV +
%              nodal GTVs
% 8. nBatch: Number of parallel batch.
%            --> Ex: 8
% 9. matlabPATH: Full path to the MATLAB excutable on the system.
%                --> Ex: 'matlab' (symbolic link)
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

startpath = pwd;
ind = strfind(namePT{1},'_'); cohortID = namePT{1}(1:ind-1);
cd(pathNonText), mkdir(['batchLog_',cohortID,'_',featType]), cd(['batchLog_',cohortID,'_',featType]), pathBatch = pwd;
time = 10; % Number of seconds to wait before checking if parallel computations are done
nameOutcomes = fieldnames(outcomes); nOutcomes = numel(nameOutcomes);

features = {'SUVmax','SUVpeak','SUVmean','aucCSH','TLG','gETU','Volume','Size','Solidity','Eccentricity'};
nFeatures = numel(features);

% PRODUCE BATCH COMPUTATIONS
nPatient = numel(namePT); valid = ones(nPatient,1);
for i = 1:nPatient
    if isempty(nameROI{i})
        valid(i) = 0;
    end
end
indValid = find(valid); nPatient = numel(indValid);
if nPatient < nBatch
    nBatch = nPatient;
end
[patients] = batchPatients(nPatient,nBatch);
save('workspace','pathData','pathNonText','namePT','nameCT','nameROI','patients','indValid','featType'), pause(5);
for i = 1:nBatch
    nameScript = ['batch',num2str(i),'_script.m'];
    fid = fopen(nameScript,'w');
    fprintf(fid,'load(''workspace'')\n');
    fprintf(fid,['calcAllNonTextureFeatures_HN_TestWithoutPercentInactive(pathData,pathNonText,namePT(indValid(patients{',num2str(i),'})),nameCT(indValid(patients{',num2str(i),'})),nameROI(indValid(patients{',num2str(i),'})),featType)\n']);
    fprintf(fid,['system(''touch batch',num2str(i),'_end'');\n']);
    fprintf(fid,'clear all');
    fclose(fid);
    system([matlabPATH,' -nojvm -nodisplay -nodesktop -nosplash < ',nameScript,' >& ',nameScript(1:end-1),'log &']);
end

% WAITING LOOP
waitBatch(pathBatch,time,nBatch)
delete('workspace.mat')

% GROUPING RESULTS FROM ALL BATCH
cd(pathNonText)
nPatient = numel(namePT);
for i = 1:nFeatures
    nonText.(features{i}).Data = zeros(nPatient,1);
    for j = 1:nPatient
        ind = strfind(namePT{j},'_'); namePatient = namePT{j}(1:ind(2)-1);
        if exist([namePatient,'_',featType,'_nonText.mat'],'file')
            tempText = load([namePatient,'_',featType,'_nonText']); tempText = struct2cell(tempText); tempText = tempText{1};
        elseif exist([namePatient,'_','GTVp','_nonText.mat'],'file')
            tempText = load([namePatient,'_','GTVp','_nonText']); tempText = struct2cell(tempText); tempText = tempText{1};
        end
        nonText.(features{i}).Data(j) = tempText.(features{i});
    end
    for j = 1:nOutcomes
        [nonText.(features{i}).Spearman.(nameOutcomes{j}).rs,nonText.(features{i}).Spearman.(nameOutcomes{j}).p] = corr(nonText.(features{i}).Data,outcomes.(nameOutcomes{j}),'type','Spearman');
    end
end
cd ..
save(['nonText_',cohortID,'_',featType],'nonText')

cd(startpath)
end