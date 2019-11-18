function calcAllAertsFeatures_batchHN_WithoutPT(pathData,pathText,nameCT,nameROI,outcomes,featType,nBatch,matlabPATH)
% -------------------------------------------------------------------------
% function calcAllSeparateTextures_batchHN(pathData,pathText,namePT,nameCT,nameROI,outcomes,featType,scale_mat,Ng_mat,nBatch,matlabPATH)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function computes SEPARATE texture features for all patients, for
% all different combinations of the following texture extraction parameters:
% - Scale: Resolution at which the ROI is isotropically resampled.
% - Ng: Number of gray-levels in the quantization process. 
%
% Different extraction parameters are passed as arrays or cells in the
% function in order to test all possible combinations. This function is 
% used for SEPARATE scans specifically. See Ref. [1,2] and 'prepareVolume.m' 
% for more details.
%
% Texture features are computed for all head and neck (HN) DICOM imaging 
% data downloaded from The Cancer Imaging Archive (TCIA) website at: 
% <http://dx.doi.org/xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx>, and first organized 
% in a 'DATA' directory using the function readAllDICOM_HN.m.  Results are 
% then saved in a folder 'TEXTURES' in the HN WORKSPACE.
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
% 1. pathData: Full path to the HN sData files directory.
%              --> Ex: '/myProject/WORKSPACE/DATA'
% 2. pathText: Full path to the HN non texture features directory.
%              --> Ex: '/myProject/WORKSPACE/FEATURES/TEXTURES'
% 3. namePT: Cell of strings of all CT sData files to read
%            --> Ex: {'HGJ_001_CT.CTscan.mat';'HGJ_022_CT.CTscan.mat'}
% 4. nameROI: Cell of strings specifying the ROI names to analyze for the
%             patients defined by "namePT" and "nameCT"
%             --> Ex: {'GTV';'GTV-P'}
% 5. outcomes: Structure specifying the status (1 or 0) for different
%              outcomes in HN cancer. Contains: outcomes.Failure, 
%              outcomes.Locoregional, outcomes.Distant. See ref.[1] for 
%              more details.
% 6. featType: Either 'GTVp' for primary GTV, or 'GTVtot' for primaty GTV +
%              nodal GTVs
%              --> Ex: 'GTVp'
% 7. scale_mat: Array vector specifying the different 'Scale' values to test.
%               --> Ex: [1,2,3,4,5]
% 8. Ng_mat: Array vector specifying the different 'Ng' values to test.
%            --> Ex: [8,16,32,64]
% 9. nBatch: Number of parallel batch.
%             --> Ex: 8
% 10.  matlabPATH: Full path to the MATLAB excutable on the system.
%      --> Ex: 'matlab' (symbolic link)
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
ind = strfind(nameCT{1},'_'); cohortID = nameCT{1}(1:ind-1);
cd(pathText), mkdir(['batchLog_',cohortID,'_',featType,'_AertsFeatures']), cd(['batchLog_',cohortID,'_',featType,'_AertsFeatures']), pathBatch = pwd;
time = 60; % Number of seconds to wait before checking if parallel computations are done
outcomes = rmfield(outcomes,{'Locoregional','Distant'}); % Radiomics signature only predicts Survival
nameOutcomes = fieldnames(outcomes); nOutcomes = numel(nameOutcomes);
scans = {'CT'}; nScans = numel(scans);

% PRODUCE BATCH COMPUTATIONS
nPatient = numel(nameCT); valid = ones(nPatient,1);
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
save('workspace','pathData','pathText','nameCT','nameROI','patients','indValid','featType'), pause(5);
for i = 1:nBatch
    nameScript = ['batch',num2str(i),'_script.m'];
    fid = fopen(nameScript,'w');
    fprintf(fid,'load(''workspace'')\n');
    fprintf(fid,['calcAllAertsFeatures_HN_WithoutPT(pathData,pathText,nameCT(indValid(patients{',num2str(i),'})),nameROI(indValid(patients{',num2str(i),'})),featType)\n']);
    fprintf(fid,['system(''touch batch',num2str(i),'_end'');\n']);
    fprintf(fid,'clear all');
    fclose(fid);
    system([matlabPATH,' -nojvm -nodisplay -nodesktop -nosplash < ',nameScript,' >& ',nameScript(1:end-1),'log &']);
end

% WAITING LOOP
waitBatch(pathBatch,time,nBatch)
delete('workspace.mat')

% GROUPING RESULTS FROM ALL BATCH
nPatient = numel(nameCT);
names = {nameCT};
for scan = 1:nScans
    cd(pathText)
    sign = struct;
    signTemp = zeros(nPatient,4); % 4-feature signature
    tempText = cell(1,nPatient); % Cell used to load patient textures only once
    for p = 1:nPatient
        ind = strfind(names{scan}{p},'.'); namePatientScan = names{scan}{p}(1:ind(1)-1);
        if exist([namePatientScan,'_',featType,'_sign.mat'],'file')
            load([namePatientScan,'_',featType,'_sign']) % Variable 'signature' is now in MATLAB workspace
        else
            load([namePatientScan,'_','GTVp','_sign']) % Variable 'signature' is now in MATLAB workspace
        end
        tempText{p} = signature;
    end
    for p = 1:nPatient
        signTemp(p,:) = tempText{p};
    end
    sign.Data = signTemp;
    for o = 1:nOutcomes
        [sign.Spearman.(nameOutcomes{o}).rs,sign.Spearman.(nameOutcomes{o}).p] = corr(sign.Data,outcomes.(nameOutcomes{o}),'type','Spearman','rows','pairwise');
    end
    cd .., save(['AertsSign_',cohortID,'_',scans{scan},'_',featType],'sign')
end

cd(startpath)
end