function calcAllSeparateTextures_batchHN(pathData,pathText,namePT,nameCT,nameROI,outcomes,featType,scale_mat,algo_cell,Ng_mat,nBatch,matlabPATH)
% -------------------------------------------------------------------------
% function calcAllSeparateTextures_batchHN(pathData,pathText,namePT,nameCT,nameROI,outcomes,featType,scale_mat,Ng_mat,nBatch,matlabPATH)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function computes SEPARATE texture features for all patients, for
% all different combinations of the following texture extraction parameters:
% - Scale: Resolution at which the ROI is isotropically resampled.
% - Quantization algorith: Type of quantization algorithm used.
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
%              --> Ex: 'GTVp'
% 8. scale_mat: Array vector specifying the different 'Scale' values to test.
%               --> Ex: [1,2,3,4,5]
% 9. Ng_mat: Array vector specifying the different 'Ng' values to test.
%            --> Ex: [8,16,32,64]
% 10. nBatch: Number of parallel batch.
%             --> Ex: 8
% 11.  matlabPATH: Full path to the MATLAB excutable on the system.
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
ind = strfind(namePT{1},'_'); cohortID = namePT{1}(1:ind-1);
cd(pathText), mkdir(['batchLog_',cohortID,'_',featType,'_SepText']), cd(['batchLog_',cohortID,'_',featType,'_SepText']), pathBatch = pwd;
time = 60; % Number of seconds to wait before checking if parallel computations are done
nameOutcomes = fieldnames(outcomes); nOutcomes = numel(nameOutcomes);
scans = {'PT','CT'}; nScans = numel(scans);

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
save('workspace','pathData','pathText','namePT','nameCT','nameROI','patients','indValid','featType','scale_mat','algo_cell','Ng_mat'), pause(5);
for i = 1:nBatch
    nameScript = ['batch',num2str(i),'_script.m'];
    fid = fopen(nameScript,'w');
    fprintf(fid,'load(''workspace'')\n');
    fprintf(fid,['calcAllSeparateTextures_HN(pathData,pathText,namePT(indValid(patients{',num2str(i),'})),nameCT(indValid(patients{',num2str(i),'})),nameROI(indValid(patients{',num2str(i),'})),featType,scale_mat,algo_cell,Ng_mat)\n']);
    fprintf(fid,['system(''touch batch',num2str(i),'_end'');\n']);
    fprintf(fid,'clear all');
    fclose(fid);
    system([matlabPATH,' -nojvm -nodisplay -nodesktop -nosplash < ',nameScript,' >& ',nameScript(1:end-1),'log &']);
end

% WAITING LOOP
waitBatch(pathBatch,time,nBatch)
delete('workspace.mat')

% GROUPING RESULTS FROM ALL BATCH
nPatient = numel(namePT);
names = {namePT,nameCT};
for scan = 1:nScans
    cd(pathText)
    if exist(['HGJ_001_',scans{scan},'_',featType,'_text.mat'],'file')
        temp = load(['HGJ_001_',scans{scan},'_',featType,'_text']); temp = struct2cell(temp); temp = temp{1}; % In order to get the necessary 'nameType' and 'nameFeature' fields
    else
        temp = load(['HGJ_001_',scans{scan},'_','GTVp','_text']); temp = struct2cell(temp); temp = temp{1}; % In order to get the necessary 'nameType' and 'nameFeature' fields
    end
    nameType = fieldnames(temp.Experiment1); nameType(end) = []; nType = numel(nameType); % All texture types are the same
    text = cell(numel(scale_mat),numel(algo_cell),numel(Ng_mat));
    tempText = cell(1,nPatient); % Cell used to load patient textures only once
    for p = 1:nPatient
        ind = strfind(names{scan}{p},'.'); namePatientScan = names{scan}{p}(1:ind(1)-1);
        if exist([namePatientScan,'_',featType,'_text.mat'],'file')
            load([namePatientScan,'_',featType,'_text']) % Variable 'textures' is now in MATLAB workspace
        else
            load([namePatientScan,'_','GTVp','_text']) % Variable 'textures' is now in MATLAB workspace
        end
        tempText{p} = textures;
    end
    experiment = 0;
    for s = 1:numel(scale_mat)
        for a = 1:numel(algo_cell)
            for n = 1:numel(Ng_mat)
                text{s,a,n} = struct;
                experiment = experiment + 1;
                strExperiment = ['Experiment',num2str(experiment)];
                for t = 1:nType
                    nameFeature = fieldnames(temp.(strExperiment).(nameType{t})); nFeature = numel(nameFeature);
                    for f = 1:nFeature
                        data = zeros(nPatient,1);
                        for p = 1:nPatient
                            data(p,1) = tempText{p}.(strExperiment).(nameType{t}).(nameFeature{f});
                        end
                        text{s,a,n}.(nameType{t}).(nameFeature{f}).Data = data;
                        for o = 1:nOutcomes
                            [text{s,a,n}.(nameType{t}).(nameFeature{f}).Spearman.(nameOutcomes{o}).rs,text{s,a,n}.(nameType{t}).(nameFeature{f}).Spearman.(nameOutcomes{o}).p] = corr(data,outcomes.(nameOutcomes{o}),'type','Spearman');
                        end
                    end
                end
            end
        end
    end
    cd .., save(['text_',cohortID,'_',scans{scan},'_',featType],'text')
end

cd(startpath)
end