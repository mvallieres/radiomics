function readAllTCGA_LGG(pathDICOM,pathROI,pathIMAGES,patientID,scans,T1Wpath,T1CEpath,T2Wpath,T2Fpath)
% -------------------------------------------------------------------------
% function readAllTCGA_LGG(pathDICOM,pathROI,pathIMAGES,patientID,T1Wpath,T1CEpath,T2Wpath,T2Fpath)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% Read all DICOM images and ROI masks downloaded from the TCGA website at:
% http://doi.org/10.7937/K9/TCIA.2016.L4LTD3TK. This function converts the
% DICOM images and ROI masks into 'key' variables for the available 'T1W', 
% 'T1CE', 'T2W' and 'T2F' imaging series of the TCGA patients defined in 
% the 'patientID' input. See 'OUTPUTS' below for a description of the 'key'
% variable.
% -------------------------------------------------------------------------
% INPUTS:
% 1. pathDICOM: Full path to the 'DOI' folder containing the DICOM images
%               of the LGG data from the TCIA website.
%               --> Ex: '/myProject/WORKSPACE/TCGA_DATA/DOI'
% 2. pathROI: Full path to the 'ROI' folder containing the DICOM images
%             of the LGG data from the TCIA website
%             --> Ex: '/myProject/WORKSPACE/TCGA_DATA/ROI'
% 3. pathIMAGES: Full path to the 'IMAGES' folder where we will save the
%                'key' variables for each imaging series containing the 
%                DICOM images and ROI definition in MATLAB format.
%                --> Ex: '/myProject/WORKSPACE/TCGA_DATA/IMAGES'
% 4. patientID: Cell of strings of size {nPatient X 1} specifying the IDs
%               of the LGG patients used in this study.
%               --> Ex: {108 X 1} cell of strings
% 5. scans: Cell of strings specifying the type of MRI scans used in this
%           study.
%           --> Ex: {'T1W','T1CE','T2W','T2F'}
% 6. T1Wpath: Cell of strings of size {nPatient X 1} specifying the paths
%             (SeriesInstanceUID) of the T1-weighted images used in this
%             study.
%             --> Ex: {108 X 1} cell of strings
% 7. T1CEpath: Cell of strings of size {nPatient X 1} specifying the paths
%              (SeriesInstanceUID) of the T1-weighted post-contrast images 
%              used in this study.
%              --> Ex: {108 X 1} cell of strings
% 8. T2Wpath: Cell of strings of size {nPatient X 1} specifying the paths
%             (SeriesInstanceUID) of the T2-weighted images used in this
%             study.
%             --> Ex: {108 X 1} cell of strings
% 9. T2Fpath: Cell of strings of size {nPatient X 1} specifying the paths
%             (SeriesInstanceUID) of the T2-flair images used in this
%             study.
%             --> Ex: {108 X 1} cell of strings
% -------------------------------------------------------------------------
% OUTPUTS: Multiple 'key' variable files for each imaging series of each
%          patient saved in 'pathIMAGES'.
%          --> key.volume: 3D array specifying the imaging volume
%          --> key.mask: 3D array specifying the ROI mask in key.volume
%          --> key.dims: Structure specifying geometrical dimensions of the
%              voxels in key.volume and key.mask
% -------------------------------------------------------------------------
% AUTHOR(S): 
% - Martin Vallieres <mart.vallieres@gmail.com>
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

%%%% IMPORTANT!!! %%% 
% ERROR NOTICED WHEN CREATING DICOM ROI MASKS
% When reading DICOM data from TCIA:
% --> Don't read image 000048.dcm of "TCGA-FG-7643", SERIES
% "1.3.6.1.4.1.14519.5.2.1.5826.4003.310945537854395744261987339837".
% THAT IMAGE IS NOT THE SAME SIZE AS THE OTHER IMAGES.
% So a temporary specialized function dedicated to reading the DICOM data
% will be create for this imaging series of this patient.

end