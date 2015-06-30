function calcAllFeatureSets_STS(pathWORK,pathMINE,fSetNameType,outcome,setSize,nonTextStruct,textCells,textCellsName,paramAll,freedomMat,baseline,alpha,delta,nBoot)
% -------------------------------------------------------------------------
% function calcAllFeatureSets_STS(pathWORK,pathMINE,fSetNameType,outcome,setSize,nonTextStruct,textCells,textCellsName,paramAll,freedomMat,baseline,alpha,delta,nBoot)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function computes feature set reduction for a given feature set type 
% and for all experiments with different degrees of freedom. See ref. [1] 
% for more details.
% -------------------------------------------------------------------------
% REFERENCE:
% [1] Vallieres, M. et al. (2015). A radiomics model from joint FDG-PET and 
%     MRI texture features for the prediction of lung metastases in soft-tissue 
%     sarcomas of the extremities. Physics in Medicine and Biology, 60(14), 
%     5471-5496. doi:10.1088/0031-9155/60/14/5471
% -------------------------------------------------------------------------
% INPUTS:
% - pathWORK: Full path to the STS WORKSPACE directory.
% - pathMINE: Full path to the MINE.jar executable. The executable can be
%             downloaded at: <http://www.exploredata.net/>.
% - fSetNameType: String specifying the name of the type of feature set 
%                 (e.g., 'PET', 'SEPARATE', 'FUSED', etc.)
% - outcome: Column vector of size [nInst X 1] specifying the outcome status 
%            (1 or 0) for all instances.
% - setSize: Size of the output feature set (typically set to 25 in ref. [1]).
% - nonTextStruct: Structure data for non-texture features. This structure 
%                  is of the same format as the one saved as output to 
%                  computeAllNonTextureFeatures_STS.m, for example.                                
% - textCells: Cell vector of organized texture cells for all the
%              different scans. Format: textCells = {cell1, cell2, etc.}, 
%              where cell1, cell2, etc. are the files saved as output to
%              organizeSeparateTextures_STS.m or organizeFusedTextures_STS.m,
%              for example.
% - textCellsName: Cell vector corresponding to the name of the
%                  corresponding cells in textCells. 
%                  (e.g., textCellsName = {'PET','T1','T2FS'} for SEPARATE
%                  scans,  textCellsName = {'PET_T1','PET_T2FS'} for FUSED 
%                  scans, as defined in ref. [1])
% - paramAll: Cell vector incorporating all texture extraction parameters
%             tested in textCells. See EXAMPLE below for more details.
% - freedomMat:  Matrix of row vectors of 1's and 0's to specify the degree 
%                of freedom on texture extraction parameters for all 
%                experiments. For example, for an nth experiment where 
%                extraction parameters 1, 2 and 4 in paramAll are allowed 
%                to vary, use freedomMat(n,:) = [1,1,0,1].
% - baseline: Vector of numerical values specifying the baseline texture 
%             extraction parameters for each entry in paramAll. See EXAMPLE
%             below for more details.
% - alpha: Numerical values specifying the coefficient of the first part of
%          the Gain equation, as defined in ref. [1].
% - delta: Numerical values specifying the coefficient of the second part 
%          of the Gain equation, as defined in ref. [1] (third part is set
%          to 0 in this function).
% - nBoot: Number of bootstrap samples to use.
%
% See <https://github.com/mvallieres/radiomics/tree/master/STUDIES/STSstudy/Functions>
% to find computeAllNonTextureFeatures_STS.m, organizeSeparateTextures_STS.m
% and organizeFusedTextures_STS.m. See masterScript_STS.m for a complete 
% example of how to utilize the current function.
% -------------------------------------------------------------------------
% OUTPUTS: Feature sets are saved in a folder named 'FSET' in the STS
%          WORKSPACE.
% -------------------------------------------------------------------------
% EXAMPLE:
% MRIinv_cell = {'NoInv','Inv'};
% MRIweight_mat = [1/4,1/3,1/2,2/3,3/4];
% R_mat = [1/2,2/3,1,3/2,2];
% scale_cell = {'pixelW',1,2,3,4,5};
% algo_cell = {'Equal','Lloyd'};
% Ng_mat = [8,16,32,64];
%
% FOR FUSED SCANS
% paramAll = {MRIinv_cell,MRIweight_mat,R_mat,scale_cell,algo_cell,Ng_mat};
% freedomMat = [1 1 1 1 ; 0 0 0 0 ; 0 1 0 1]; (example)
% baseline = [1 3 3 1 2 3];
%
% FOR SEPARATE SCANS
% paramAll = {R_mat,scale_cell,algo_cell,Ng_mat};
% freedomMat = [1 1 1 1 1 1 ; 0 0 0 0 0 0 ;  0 1 0 1 0 1]; (example)
% baseline = [3 1 2 3];
%
% NOTE: paramAll must always be of the same format and size for SEPARATE 
% and FUSED SCANS, with the same ordering of different extraction
% parameters.
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

startpath = pwd;
cd([pathWORK,'/FSET'])

nParamType = size(freedomMat,2);
nFreedom = size(freedomMat,1);
tStart = tic;
for i = 1:nFreedom
    nameSave=['FSET_',fSetNameType,'_'];
    for j = 1:nParamType
        nameSave = [nameSave,num2str(freedomMat(i,j))];
    end
    fprintf(['COMPUTING ',nameSave,' ... '])
    tic
    [fSet] = featureSetReduction(pathMINE,fSetNameType,outcome,setSize,nonTextStruct,textCells,textCellsName,paramAll,freedomMat(i,:),baseline,alpha,delta,nBoot);
    toc
    save(nameSave,'fSet')
end
time = toc(tStart);
fprintf('TOTAL TIME: %.2f seconds\n',time)

cd(startpath)
end