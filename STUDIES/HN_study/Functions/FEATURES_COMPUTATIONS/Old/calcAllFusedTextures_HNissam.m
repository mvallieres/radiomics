function calcAllFusedTextures_HNissam(pathWORK,nPatient,roiNumb,CTinv_cell,CTweight_mat,R_mat,scale_cell,algo_cell,Ng_mat,fusWeights)
% -------------------------------------------------------------------------
% function calcAllFusedTextures_HNissam(pathWORK,nPatient,roiNumb,CTinv_cell,CTweight_mat,R_mat,scale_cell,algo_cell,Ng_mat)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function computes FUSED texture features for all patients, for all 
% different combinations of the following texture extraction parameters:
% - CT Inv.: Inversion of CT intensities in the PET/CT fusion process.
% - CT weight: Weight given to MRI wavelet low-pass sub-bands in the 
%               PET/CT fusion process. 
% - WPBF ratio 'R': Ratio specifying the ratio of weight given to band-pass 
%                   coefficients to the weigth given to the rest of 
%                   coefficients (HHH and LLL) in the wavelet band-pass 
%                   filtering process.
% - Scale: Resolution at which the ROI is isotropically resampled.
% - Quant. Algo.: Quantization algorithm used on analyzed ROI.
% - Ng: Number of gray-levels in the quantization process. 
%
% Different extraction parameters are passed as arrays or cells in the
% function in order to test all possible combinations. This function is 
% used for FUSED scans specifically. See Ref. [1,2] and 'prepareVolume.m' 
% for more details.
%
% Texture features are computed for all head and neck (HN) DICOM  imaging 
% data downloaded from The Cancer Imaging Archive (TCIA) website at: 
% <http://dx.doi.org/10.7937/K9/xxxxxxxxxxxxxxxxxxx, and first organized 
% in a 'DATA' directory using the function readAllDICOM_HN.m.  Results are 
% then saved in a folder 'TEXTURES' in the HN WORKSPACE.
% -------------------------------------------------------------------------
% REFERENCE:
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
% - nPatient: Number of patients to read.
% - roiNumb: Array of size [N X 2] specifying the contour number 
%            corresponding to the GTV in sData{2}.scan.contour for all 
%            patients. First column is for the PET data, second column for 
%            CT data.
% - CTinv_cell: Cell vector specifying the different 'CT Inv.' values to test.
% - CTweight_mat: Array vector specifying the different 'CT weight' values to test.
% - R_mat: Array vector specifying the different 'R' ratio values to test.
% - scale_cell: Cell vector specifying the different 'Scale' values to test.
% - algo_cell: Cell vector specifying the different 'Quant. Algo.' values to test.
% - Ng_mat: Array vector specifying the different 'Ng' values to test.
% -------------------------------------------------------------------------
% EXAMPLE:
% CTinv_cell = {'NoInv','Inv'};
% CTweight_mat = [1/4,1/3,1/2,2/3,3/4];
% R_mat = [1/2,2/3,1,3/2,2];
% scale_cell = {'pixelW',1,2,3,4,5};
% algo_cell = {'Uniform','Equal','Lloyd'};
% Ng_mat = [8,16,32,64];
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

startpath=pwd;
cd(pathWORK), cd('TEXTURES'), pathTEXT = pwd;
cd .., cd('DATA'), pathDATA = pwd;

% INITIALIZATION
scans = {'CT'};
nScans = numel(scans);
suffix = ['_TEXT'];

% COMPUTATION
fprintf('\n')
for i = 1:nPatient
    load(['Patient',num2str(i),'_PET']), sDataPET = sData;
    for j = 1:nScans
        nameLoad = ['Patient',num2str(i),'_CT'];
        load(nameLoad), sDataCT = sData;
        tStart = tic;
        fprintf(['\n*********************** COMPUTING TEXTURES: PATIENT %u, FUSED PET/',scans{j},' SCAN ***********************'],i)
        [textures] = calcPatientFusText_HNissam(sDataPET,sDataCT,roiNumb(i,:),CTinv_cell,CTweight_mat,R_mat,scale_cell,algo_cell,Ng_mat,fusWeights);
        cd(pathTEXT), save(['Patient',num2str(i),'_PET_',scans{j},suffix],'textures')
        cd(pathDATA)
        time = toc(tStart);
        fprintf('TOTAL TIME: %.2f seconds\n',time)
    end
end

cd(startpath)
end