function calcAllAertsFeatures_HN(pathData,pathText,namePT,nameCT,nameROI,featType)
% -------------------------------------------------------------------------
% function calcAllSeparateTextures_HN(pathData,pathText,namePT,nameCT,nameROI,featType,scale_mat,Ng_mat)
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
% 2. pathNonText: Full path to the HN non texture features directory.
%              --> Ex: '/myProject/WORKSPACE/FEATURES/NON_TEXTURES'
% 3. namePT: Cell of strings of all PET sData files to read
%            --> Ex: {'HGJ_001_PT.PTscan.mat';'HGJ_022_PT.PTscan.mat'}
% 4. namePT: Cell of strings of all CT sData files to read
%            --> Ex: {'HGJ_001_CT.CTscan.mat';'HGJ_022_CT.CTscan.mat'}
% 5. nameROI: Cell of strings specifying the ROI names to analyze for the
%             patients defined by "namePT" and "nameCT"
%             --> Ex: {'GTV';'GTV-P'}
% 6. featType: Either 'GTVp' for primary GTV, or 'GTVtot' for primary GTV +
%              nodal GTVs
% 7. scale_mat: Array vector specifying the different 'Scale' values to test.
%               --> Ex: [1,2,3,4,5]
% 8. Ng_mat: Array vector specifying the different 'Ng' values to test.
%            --> Ex: [8,16,32,64]
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

startpath=pwd;

% INITIALIZATION
scans = {'PT','CT','CTorig'}; nScans = numel(scans); names = {namePT,nameCT,nameCT};
scanTypes = {'PET','CT','CT'};
nPatient = numel(namePT);

% COMPUTATION
fprintf('\n')
for i = 1:nPatient
    for j = 1:nScans
        cd(pathData)
        tStart = tic;
        if strcmp(scans{j},'CTorig')
            fprintf(['\n*********************** COMPUTING AERTS SIGNATURE: %s, ORIGINAL ***********************'],names{j}{i});
        else
            fprintf(['\n*********************** COMPUTING AERTS SIGNATURE: %s, ***********************'],names{j}{i});
        end
        scanType = scanTypes{j};
        load(names{j}{i}) % Variable 'sData' now in MATLAB Workspace;
        spatialRef = sData{2}.scan.volume.spatialRef;
        mask = zeros(size(sData{2}.scan.volume.data));
        [contourVect] = findContours_HN(sData,{nameROI{i}},scanType);
        for c = 1:numel(contourVect)
            temp = zeros(size(sData{2}.scan.volume.data));
            [~,ROImask,bb] = getROI(sData,contourVect(c));
            temp(bb(1,1):bb(1,2),bb(2,1):bb(2,2),bb(3,1):bb(3,2)) = ROImask;
            mask = mask + temp;
        end
        mask(mask > 0) = 1;
        [bb] = computeBoundingBox(mask); a = 0; b = 0; c = 0;
        sizeROI = [bb(1,2)-bb(1,1)+1,bb(2,2)-bb(2,1)+1,bb(3,2)-bb(3,1)+1];
        if mod(sizeROI(1),2)
            a = a + 1;
        end
        if mod(sizeROI(2),2)
            b = b + 1;
        end
        if mod(sizeROI(3),2)
            c = c + 1;
        end
        mask = mask(bb(1,1):bb(1,2)+a,bb(2,1):bb(2,2)+b,bb(3,1):bb(3,2)+c);
        vol = sData{2}.scan.volume.data(bb(1,1):bb(1,2)+a,bb(2,1):bb(2,2)+b,bb(3,1):bb(3,2)+c);
        if strcmp(scanType,'PET')
            [vol] = computeSUVmap(vol,sData{3}(1)); 
        end
        if strcmp(scans{j},'CTorig')
            [signature] = getAertsSign(vol,mask,spatialRef,scanType); % Calculating the radiomic signature as originially defined by Aerts et al.
            suffix = '_orig';
        else
            [signature] = getAertsSignMod(vol,mask,spatialRef,scanType); % Improved version (dimensionless compactness, 1 GLRLM matrix, GLN normalization)
            suffix = '';
        end
        ind = strfind(names{j}{i},'.'); namePatientScan = names{j}{i}(1:ind(1)-1);
        cd(pathText), save([namePatientScan,'_',featType,'_sign',suffix],'signature')
        time = toc(tStart);
        fprintf('TOTAL TIME: %.2f seconds\n',time)
    end
end

cd(startpath)
end