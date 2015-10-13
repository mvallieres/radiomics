function [aucCSH] = getAUCCSH(ROIonlyPET)
% -------------------------------------------------------------------------
% function [aucCSH] = getAUCCSH(ROIonlyPET)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function computes the area under the curve of the cumulative 
% SUV-volume histogram from the region of interest (ROI) of an input PET 
% volume.
% -------------------------------------------------------------------------
% REFERENCE:
% [1] van Velden, F. H. P. et al. (2011). Evaluation of a cumulative
%     SUV-volume histogram method for parameterizing heterogeneoous
%     intratumoural FDG uptake in non-small cell lung cancer. European 
%     Journal of Nuclear Medicine and Molecular Imaging, 38(9), 1636-1647.
% -------------------------------------------------------------------------
% INPUTS:
% - ROIonlyPET: 3D array representing the PET volume in SUV format, with 
%               voxels outside the ROI set to NaNs. 
% -------------------------------------------------------------------------
% OUTPUTS:
% - aucCSH: Area under the curve of the cumulative SUV-volume histogram
%           describing the percentage of total volume of the ROI above a 
%           percentage threshold of maximum SUV.
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

nBins = 1000; % By default.

outliers = find(ROIonlyPET > mean(ROIonlyPET(~isnan(ROIonlyPET(:)))) + 3*std(ROIonlyPET(~isnan(ROIonlyPET(:)))));
goodVoxels = find(ROIonlyPET <= mean(ROIonlyPET(~isnan(ROIonlyPET(:)))) + 3*std(ROIonlyPET(~isnan(ROIonlyPET(:)))));
ROIonlyPET(outliers) = mean(ROIonlyPET(goodVoxels));

ROIonlyPET = ROIonlyPET-min(ROIonlyPET(:));
ROIonlyPET = ROIonlyPET./max(ROIonlyPET(:));
volume = ROIonlyPET(~isnan(ROIonlyPET));
nVoxel = numel(volume);

bins = hist(volume,nBins);
CSH = fliplr(cumsum(fliplr(bins))./nVoxel);
aucCSH = sum(CSH./nBins);

end