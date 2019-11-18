function calcAllNonTextureFeatures_HN_TestWithoutPercentInactive(pathData,pathNonText,namePT,nameCT,nameROI,featType)
% -------------------------------------------------------------------------
% function calcAllNonTextureFeatures_HN(pathData,pathNonText,namePT,nameCT,nameROI)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% Computes 11 non-texture features for all head and neck (HN) DICOM 
% imaging data downloaded from The Cancer Imaging Archive (TCIA) website 
% at: <http://dx.doi.org/xxxxxxxxxxxxxxxxxxxxxxxxxxxxx>, and organized in a 
% 'DATA' directory using the function readAllDICOM_HN.m.
% -------------------------------------------------------------------------
% REFERENCE:
% [1] Vallieres, M. et al. (2016). FDG-PET/CT radiomics models for the 
%     early prediction of different tumour outcomes in head and neck cancer.
%     The Journal of Nuclear Medicine, aa(bb), xxx-yyy. 
%     doi:
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
% 6. featType: Either 'GTVp' for primary GTV, or 'GTVtot' for primaty GTV +
%              nodal GTVs
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

% INTIALIZATION
nonTextures = struct;
nPatient = numel(namePT);


% COMPUTATION OF FEATURES
fprintf('\n')
for i = 1:nPatient
    cd(pathData)
    
    % SUV metrics
    load(namePT{i}) % Variable 'sData' now in MATLAB Workspace
    mask = zeros(size(sData{2}.scan.volume.data));
    [contourVect] = findContours_HN(sData,{nameROI{i}});
    for c = 1:numel(contourVect)
        temp = zeros(size(sData{2}.scan.volume.data));
        [~,ROImask,bb] = getROI(sData,contourVect(c));
        temp(bb(1,1):bb(1,2),bb(2,1):bb(2,2),bb(3,1):bb(3,2)) = ROImask;
        mask = mask + temp;
    end
    mask(mask > 0) = 1; [bb] = computeBoundingBox(mask);
    SUVmap = sData{2}.scan.volume.data(bb(1,1):bb(1,2),bb(2,1):bb(2,2),bb(3,1):bb(3,2));
    ROImask = mask(bb(1,1):bb(1,2),bb(2,1):bb(2,2),bb(3,1):bb(3,2));
    [SUVmap] = computeSUVmap(SUVmap,sData{3}(1)); SUVmap(~ROImask) = NaN;
    fprintf('COMPUTING ''SUV METRICS'' of %s ... ',namePT{i})
    [SUVmax,SUVpeak,SUVmean,aucCSH] = getSUVmetrics(SUVmap);
    [percentInactive] = getPercentInactive(SUVmap,0.0065);
    [Ma_gETU]= getMa_gETU(SUVmap,0.25);
    fprintf('DONE\n')
    nonTextures.SUVmax = SUVmax;
    nonTextures.SUVpeak = SUVpeak;
    nonTextures.SUVmean = SUVmean;
    nonTextures.aucCSH = aucCSH;
    %nonTextures.PercentInactive = percentInactive;
    
    % Volume, Size, Solidity, Eccentricity
    load(nameCT{i}) % Variable 'sData' now in MATLAB Workspace
    mask = zeros(size(sData{2}.scan.volume.data));
    [contourVect] = findContours_HN(sData,{nameROI{i}});
    for c = 1:numel(contourVect)
        temp = zeros(size(sData{2}.scan.volume.data));
        [~,ROImask,bb] = getROI(sData,contourVect(c));
        temp(bb(1,1):bb(1,2),bb(2,1):bb(2,2),bb(3,1):bb(3,2)) = ROImask;
        mask = mask + temp;
    end
    mask(mask > 0) = 1; [bb] = computeBoundingBox(mask);
    ROIbox = sData{2}.scan.volume.data(bb(1,1):bb(1,2),bb(2,1):bb(2,2),bb(3,1):bb(3,2));
    ROImask = mask(bb(1,1):bb(1,2),bb(2,1):bb(2,2),bb(3,1):bb(3,2));
    pixelW = sData{2}.scan.volume.spatialRef.PixelExtentInWorldX; sliceS = sData{2}.scan.volume.spatialRef.PixelExtentInWorldZ;
    ROIonly = ROIbox; ROIonly(~ROImask) = NaN;
    fprintf(['COMPUTING ''VOLUME'' of ',nameCT{i},' ... '])
    [volumeROI] = getVolume(ROIonly,pixelW,sliceS); volumeROI = volumeROI/1000; % in cm^3
    nonTextures.Volume = volumeROI;
    fprintf('DONE\n')
    fprintf(['COMPUTING ''TLG'' of ',namePT{i},' ... '])
    nonTextures.TLG = SUVmean * volumeROI;
    fprintf('DONE\n')
    fprintf(['COMPUTING ''gETU'' of ',namePT{i},' ... '])
    nonTextures.gETU = Ma_gETU * volumeROI^(1/0.25);
    fprintf('DONE\n')
    fprintf(['COMPUTING ''SIZE'' of ',nameCT{i},' ... '])
    [sizeROI] = getSize(ROIonly,pixelW,sliceS); nonTextures.Size = sizeROI;
    fprintf('DONE\n')
    fprintf(['COMPUTING ''SOLIDITY'' of ',nameCT{i},' ... '])
    [solidity] = getSolidity(ROIonly,pixelW,sliceS); nonTextures.Solidity = solidity;
    fprintf('DONE\n')
    fprintf(['COMPUTING ''ECCENTRICITY'' of ',nameCT{i},' ... '])
    [eccentricity] = getEccentricity(ROIonly,pixelW,sliceS); nonTextures.Eccentricity = eccentricity;
    fprintf('DONE\n')
    fprintf('\n')
    
    cd(pathNonText)
    ind = strfind(namePT{i},'_'); namePatient = namePT{i}(1:ind(2)-1);
    save([namePatient,'_',featType,'_nonText'],'nonTextures')
end

cd(startpath)
end