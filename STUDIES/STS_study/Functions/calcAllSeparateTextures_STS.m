function calcAllSeparateTextures_STS(pathWORK,nPatient,roiNumb,R_mat,scale_cell,algo_cell,Ng_mat)
% -------------------------------------------------------------------------
% function calcAllSeparateTextures_STS(pathWORK,nPatient,roiNumb,R_mat,scale_cell,algo_cell,Ng_mat)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function computes SEPARATE texture features for all patients, for
% all different combinations of the following texture extraction parameters:
% - WPBF ratio 'R': Ratio specifying the ratio of weight given to band-pass 
%                   coefficients to the weigth given to the rest of 
%                   coefficients (HHH and LLL) in the wavelet band-pass 
%                   filtering process.
% - Scale: Resolution at which 'volume' is isotropically resampled.
% - Quant. Algo.: Quantization algorithm used on 'volume'.
% - Ng: Number of gray-levels on the quantization process. 
%
% Different extraction parameters are passed as arrays or cells in the
% function in order to test all possible combinations. This function is 
% used for SEPARATE scans specifically. See Ref. [1] and 'prepareVolume.m' 
% for more details.
%
% Texture features are computed for all soft-tissue sarcoma (STS) DICOM 
% imaging data downloaded from The Cancer Imaging Archive (TCIA) website 
% at: <http://dx.doi.org/10.7937/K9/TCIA.2015.7GO2GSKS>, and first organized 
% in a 'DATA' directory using the function readAllDICOM_STS.m.  Results are 
% then saved in a folder 'TEXTURES' in the STS WORKSPACE.
% -------------------------------------------------------------------------
% REFERENCE:
% [1] Vallieres, M. et al. (2015). A radiomics model from joint FDG-PET and 
%     MRI texture features for the prediction of lung metastases in soft-tissue 
%     sarcomas of the extremities. Physics in Medicine and Biology, 60(14), 
%     5471-5496. doi:10.1088/0031-9155/60/14/5471
% -------------------------------------------------------------------------
% INPUTS:
% - pathWORK: Full path to the STS WORKSPACE directory.
% - nPatient: Number of patients to read.
% - roiNumb: Vector of the contour number to use in sData{2}.scan.contour
%            for all patients. (Load 'contour_Mass' from the STS WORKSPACE)
% - R_mat: Array vector specifying the different 'R' ratio values to test.
% - scale_cell: Cell vector specifying the different 'Scale' values to test.
% - algo_cell: Cell vector specifying the different 'Quant. Algo.' values to test.
% - Ng_mat: Array vector specifying the different 'Ng' values to test.
% -------------------------------------------------------------------------
% EXAMPLE:
% R_mat = [1/2,2/3,1,3/2,2];
% scale_cell = {'pixelW',1,2,3,4,5};
% algo_cell = {'Equal','Lloyd'};
% Ng_mat = [8,16,32,64];
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

startpath=pwd;
cd(pathWORK), cd('TEXTURES'), pathTEXT = pwd;
cd .., cd('DATA'), pathDATA = pwd;

% INITIALIZATION
scans = {'PET','T1','T2FS'};
nScans = numel(scans);
suffix = ['_TEXT'];

% COMPUTATION
fprintf('\n')
for i = 1:nPatient
    for j = 1:nScans
        cd(pathDATA)
        try
            nameLoad = ['Patient',num2str(i),'_',scans{j}];
            load(nameLoad) % Variable 'sData' now in MATLAB Workspace
        catch
            nameLoad = ['Patient',num2str(i),'_STIR'];
            load(nameLoad) % Variable 'sData' now in MATLAB Workspace
        end
        tStart = tic;
        fprintf(['\n*********************** COMPUTING TEXTURES: PATIENT %u, ',scans{j},' SCAN ***********************'],i)
        [textures] = calcPatientSepText_STS(sData,roiNumb(i),R_mat,scale_cell,algo_cell,Ng_mat);
        cd(pathTEXT), save(['Patient',num2str(i),'_',scans{j},suffix],'textures')
        time = toc(tStart);
        fprintf('TOTAL TIME: %.2f seconds\n',time)
    end
end

cd(startpath)
end