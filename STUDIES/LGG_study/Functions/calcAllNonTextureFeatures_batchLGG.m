function calcAllNonTextureFeatures_batchLGG(pathIMAGES,pathNonText,outcomesTCGA,nBatch,matlabPATH)
% -------------------------------------------------------------------------
% function calcAllNonTextureFeatures_batchLGG(pathIMAGES,pathNonText,outcomesTCGA,nBatch,matlabPATH)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% Computes 4 non-texture features (Size, Volume, Eccentricity, Solidity for
% all 'key' variables saved in 'pathIMAGES' In addition, the 
% Spearman's rank correlation between each non-texture features and 
% different outcomes in LGG cancer. Results for each imaging series of each
% patient are saved in 'pathNonText'. Organized results are saved as a 
% structure named 'nonTextures.mat' in the TCGA_DATA folder.
% -------------------------------------------------------------------------
% INPUTS:
% 1. pathIMAGES: Full path to the 'IMAGES' folder where we will save the
%                'key' variables for each imaging series containing the 
%                DICOM images and ROI definition in MATLAB format.
%                --> Ex: '/myProject/WORKSPACE/TCGA_DATA/IMAGES'
% 2. pathNonText: Full path to the non texture feature directory where 
%                 data for each imaging volume will be saved. Defined and 
%                 created in masterScript_LGG.m
%                 --> Ex: '/myProject/WORKSPACE/TCGA_DATA/NON_TEXTURES'
% 3. outcomesTCGA: Structure specifying the status (1 or 0) for different
%                  outcomes in LGG cancer. Contains: outcomes.nonIDH1, 
%                  outcomes.IDHcodel, outcomes.progression and
%                  outcomes.lowGrade.
%                  --> structure containing 4 [nPatient X 1] vectors
% 4. nBatch: Number of parallel batch.
%            --> Ex: 8
% 5. matlabPATH: Full path to the MATLAB executable on the system.
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