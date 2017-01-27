function calcAllTextures_batchLGG(pathIMAGES,pathText,outcomesTCGA,scale_mat,algo_cell,Ng_mat,nBatch,matlabPATH)
% -------------------------------------------------------------------------
% function calcAllTextures_batchLGG(pathIMAGES,pathText,outcomesTCGA,scale_mat,algo_cell,Ng_mat,nBatch,matlabPATH)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function computes texture features for all imaging series of all 
% LGG patients, for all different combinations of the following texture 
% extraction parameters:
% - Scale: Resolution at which the ROI is isotropically resampled.
% - Quant.Algo.: Quantization algorithm used.
% - Ng: Number of gray-levels in the quantization process. 
%
% Different extraction parameters are passed as arrays or cells in the
% function in order to test all possible combinations. See Ref. [1] 
% and 'prepareVolume.m' for more details.
% -------------------------------------------------------------------------
% REFERENCES:             
% [1] Vallieres, M. et al. (2015). A radiomics model from joint FDG-PET and 
%     MRI texture features for the prediction of lung metastases in soft-tissue 
%     sarcomas of the extremities. Physics in Medicine and Biology, 60(14), 
%     5471-5496. doi:10.1088/0031-9155/60/14/5471
% -------------------------------------------------------------------------
% INPUTS:
% 1. pathIMAGES: Full path to the 'IMAGES' folder where we will save the
%                'key' variables for each imaging series containing the 
%                DICOM images and ROI definition in MATLAB format.
%                --> Ex: '/myProject/WORKSPACE/TCGA_DATA/IMAGES'
% 2. pathText: Full path to the texture feature directory where data for 
%              each imaging volume will be saved. Defined and created in 
%              masterScript_LGG.m
%              --> Ex: '/myProject/WORKSPACE/TCGA_DATA/TEXTURES'
% 3. outcomesTCGA: Structure specifying the status (1 or 0) for different
%                  outcomes in LGG cancer. Contains: outcomes.nonIDH1, 
%                  outcomes.IDHcodel, outcomes.progression and
%                  outcomes.lowGrade.
%                  --> Ex: Structure containing 4 [nPatient X 1] vectors
% 4. scale_mat: Array vector specifying the different 'Scale' values to test.
%               --> Ex: [1,2,3,4,5]
% 5. algo_cell: Cell of strings specifying the different 'Quant.Algo.' values to test.
%               --> Ex: {'Equal','Uniform'}
% 6. Ng_mat: Array vector specifying the different 'Ng' values to test.
%            --> Ex: [8,16,32,64]
% 7. nBatch: Number of parallel batch.
%             --> Ex: 8
% 8. matlabPATH: Full path to the MATLAB executable on the system.
%                --> Ex: 'matlab' (symbolic link)
% -------------------------------------------------------------------------
% AUTHOR(S): Martin Vallieres <mart.vallieres@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: January 2017
%--------------------------------------------------------------------------
% STATEMENT:
% This file is part of <https://github.com/mvallieres/radiomics/>, 
% a package providing MATLAB programming tools for radiomics analysis.
% --> Copyright (C) 2015-2017  Martin Vallieres
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


%%%%%%%%%%%%%%%%%%%%%%%
% CODE IN PREPARATION %
%%%%%%%%%%%%%%%%%%%%%%%

end