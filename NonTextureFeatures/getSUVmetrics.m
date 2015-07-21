function [SUVmax,SUVpeak,SUVmean,aucCSH] = getSUVmetrics(ROIonlyPET)
% -------------------------------------------------------------------------
% function [SUVmax,SUVpeak,SUVmean,aucCSH] = getSUVmetrics(ROIonlyPET)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function computes SUVmax, SUVpeak and SUVmean, AUC-CSH and Percent 
% Inactive metrics from the region of interest (ROI) of an input PET volume.
% -------------------------------------------------------------------------
% INPUTS:
% - ROIonlyPET: 3D array representing the PET volume in SUV format, with 
%               voxels outside the ROI set to NaNs. 
% -------------------------------------------------------------------------
% OUTPUTS:
% - SUVmax: Maximum SUV of the ROI.
% - SUVpeak: Average of the voxel with maximum SUV within the ROI and its 
%            26 connected neighbours.
% - SUVmean: Average SUV value of the ROI.
% - aucCSH: Area under the curve of the cumulative SUV-volume histogram
%           describing the percentage of total volume of the ROI above a 
%           percentage threshold of maximum SUV.
%           (van Velden et al., Eur J Nucl Med Mol Imaging 38(9), 2011).
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


% Initialization
ROIonlyPET = padarray(ROIonlyPET,[1 1 1],NaN);

% SUVmax
[SUVmax,indMax] = max(ROIonlyPET(:));

% SUVpeak (using 26 neighbors around SUVmax)
[indMaxX,indMaxY,indMaxZ] = ind2sub(size(ROIonlyPET),indMax);
connectivity = getneighbors(strel('arbitrary',conndef(3,'maximal')));
nPeak = length(connectivity);
neighborsMax = zeros(1,nPeak);
for i=1:nPeak
    neighborsMax(i) = ROIonlyPET(connectivity(i,1)+indMaxX,connectivity(i,2)+indMaxY,connectivity(i,3)+indMaxZ);
end
SUVpeak = mean(neighborsMax(~isnan(neighborsMax)));

% SUVmean
SUVmean=mean(ROIonlyPET(~isnan(ROIonlyPET)));

% AUC-CSH
[aucCSH] = getAUCCSH(ROIonlyPET);

end