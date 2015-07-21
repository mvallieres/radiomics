function calcAllNonTextureFeatures_STS(pathWORK,nPatient,roiNumb,outcome)
% -------------------------------------------------------------------------
% function calcAllNonTextureFeatures_STS(pathWORK,nPatient,roiNumb,outcome)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% Computes 9 non-texture features for all soft-tissue sarcoma (STS) DICOM 
% imaging data downloaded from The Cancer Imaging Archive (TCIA) website 
% at: <http://dx.doi.org/10.7937/K9/TCIA.2015.7GO2GSKS>, and organized in a 
% 'DATA' directory using the function readAllDICOM_STS.m. In addition, the 
% Spearman's rank correlation between each non-texture features and lung 
% metastases is computed. Results are saved as a structure named 
% 'nonTextures.mat' in the STS WORKSPACE.
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
% - outcome: Vector of lung metastases outcome status for all patients.
%            (Load 'outcome' from the STS WORKSPACE)
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
cd([pathWORK,'/DATA'])
nonText = struct;


% COMPUTATION
fprintf('\n')
for i = 1:nPatient
    
    % SUV metrics
    load(['Patient',num2str(i),'_PET']) % Variable 'sData' now in MATLAB Workspace
    ROIonly = getROIonly(sData,roiNumb(i));
    [SUVmap] = computeSUVmap(ROIonly,sData{3}(1));
    fprintf(['COMPUTING ''SUV METRICS'' of Patient',num2str(i),'...'])
    [SUVmax,SUVpeak,SUVmean,aucCSH] = getSUVmetrics(SUVmap);
    [percentInactive] = getPercentInactive(SUVmap,0.005);
    fprintf('DONE\n')
    nonText.SUVmax.Data(i,1) = SUVmax;
    nonText.SUVpeak.Data(i,1) = SUVpeak;
    nonText.SUVmean.Data(i,1) = SUVmean;
    nonText.aucCSH.Data(i,1) = aucCSH;
    nonText.PercentInactive.Data(i,1) = percentInactive;
    
    % Volume, Size, Solidity, Eccentricity
    try
        load(['Patient',num2str(i),'_T2FS']) % Variable 'sData' now in MATLAB Workspace
    catch
        load(['Patient',num2str(i),'_STIR']) % Variable 'sData' now in MATLAB Workspace
    end
    ROIonly = getROIonly(sData,roiNumb(i));
    pixelW = sData{2}.scan.pixelW; sliceS = sData{2}.scan.sliceS;
    fprintf(['COMPUTING ''VOLUME'' of Patient',num2str(i),'...'])
    [volumeROI] = getVolume(ROIonly,pixelW,sliceS); nonText.Volume.Data(i,1) = volumeROI;
    fprintf('DONE\n')
    fprintf(['COMPUTING ''SIZE'' of Patient',num2str(i),'...'])
    [sizeROI] = getSize(ROIonly,pixelW,sliceS); nonText.Size.Data(i,1) = sizeROI;
    fprintf('DONE\n')
    fprintf(['COMPUTING ''SOLIDITY'' of Patient',num2str(i),'...'])
    [solidity] = getSolidity(ROIonly,pixelW,sliceS); nonText.Solidity.Data(i,1) = solidity;
    fprintf('DONE\n')
    fprintf(['COMPUTING ''ECCENTRICITY'' of Patient',num2str(i),'...'])
    [eccentricity] = getEccentricity(ROIonly,pixelW,sliceS); nonText.Eccentricity.Data(i,1) = eccentricity;
    fprintf('DONE\n')
    fprintf('\n')
end
nonText.SUVmax.Spearman = corr(nonText.SUVmax.Data,outcome,'type','Spearman');
nonText.SUVpeak.Spearman = corr(nonText.SUVpeak.Data,outcome,'type','Spearman');
nonText.SUVmean.Spearman = corr(nonText.SUVmean.Data,outcome,'type','Spearman');
nonText.aucCSH.Spearman = corr(nonText.aucCSH.Data,outcome,'type','Spearman');
nonText.PercentInactive.Spearman = corr(nonText.PercentInactive.Data,outcome,'type','Spearman');
nonText.Volume.Spearman = corr(nonText.Volume.Data,outcome,'type','Spearman');
nonText.Size.Spearman = corr(nonText.Size.Data,outcome,'type','Spearman');
nonText.Solidity.Spearman = corr(nonText.Solidity.Data,outcome,'type','Spearman');
nonText.Eccentricity.Spearman = corr(nonText.Eccentricity.Data,outcome,'type','Spearman');


cd(pathWORK), save('nonTextures','nonText')

cd(startpath)
end