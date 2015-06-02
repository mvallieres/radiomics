% -------------------------------------------------------------------------
% DESCRIPTION: 
% This script computes parts of the experiments performed in ref. [1] 
% ('fast' version).
% -------------------------------------------------------------------------
% REFERENCE:
% [1] Vallieres, M. et al. (2015). A radiomics model from joint FDG-PET and 
%     MRI texture features for the prediction of lung metastases in soft-tissue 
%     sarcomas of the extremities. Physics in Medicine and Biology, 00(0), xxx-yyy. 
% -------------------------------------------------------------------------
% AUTHOR(S): Martin Vallieres, <mart.vallieres@gmail.com>
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


warning off

% TEXTURE EXTRACTION PARAMETERS
MRIinv_cell = {'NoInv','Inv'};
MRIweight_mat = [1/3,2/3];
R_mat = [2/3,3/2];
scale_cell = {4,5};
algo_cell = {'Equal','Lloyd'};
Ng_mat = [8,16];
nPatient = 51;
load('contour_Mass') % Variable 'roiNumb' is now in MATLAB workspace
load('outcome') % % Variable 'outcome' is now in MATLAB workspace


% 1. READ DATA DOWNLOADED FROM THE TCIA WEBSITE
readAllDICOM_STS([pwd,'/Soft-tissue Sarcoma'],51), cd ..
 

% 2. COMPUTE NON-TEXTURE FEATURES
computeAllNonTextureFeatures_STS(pwd,51,roiNumb,outcome)


% 3. COMPUTE ALL TEXTURE FEATURES
calcAllSeparateTextures_STS(pwd,51,roiNumb,R_mat,scale_cell,algo_cell,Ng_mat)
fprintf('\n')
calcAllFusedTextures_STS(pwd,51,roiNumb,MRIinv_cell,MRIweight_mat,R_mat,scale_cell,algo_cell,Ng_mat)


% 4. ORGANIZE ALL TEXTURE FEATURES
tStart = tic;
fprintf('\nORGANIZING TEXTURES FROM SEPARATE SCANS ... ')
organizeSeparateTextures_STS(pwd,51,R_mat,scale_cell,algo_cell,Ng_mat)
fprintf('DONE')
toc(tStart)

tStart = tic;
fprintf('\nORGANIZING TEXTURES FROM FUSED SCANS ... ')
organizeFusedTextures_STS(pwd,51,MRIinv_cell,MRIweight_mat,R_mat,scale_cell,algo_cell,Ng_mat)
fprintf('DONE')
toc(tStart)


% 5. PERFORM MULTIVARIABLE ANALYSIS
% IMPORTANT: Codes will be added here during the week of June 8, 2015.
% (Martin Vallières - June 1, 2015)
