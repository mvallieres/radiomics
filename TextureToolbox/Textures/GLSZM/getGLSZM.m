function [GLSZM] = getGLSZM(ROIOnly,levels)
% -------------------------------------------------------------------------
% [GLSZM] = getGLSZM(ROIOnly,levels)
% -------------------------------------------------------------------------
% DESCRIPTION:
% This function computes the Gray-Level Size Zone Matrix (GLSZM) of the 
% region of interest (ROI) of an input volume. The input volume is assumed 
% to be isotropically resampled. The zones of different sizes are computed 
% using 26-voxel connectivity.
%
% --> This function is compatible with 2D analysis (language not adapted in the text)
% -------------------------------------------------------------------------
% REFERENCE:
% [1] Thibault, G., Fertil, B., Navarro, C., Pereira, S., Cau, P., Levy, 
%     N., Mari, J.-L. (2009). Texture Indexes and Gray Level Size Zone 
%     Matrix. Application to Cell Nuclei Classification. In Pattern 
%     Recognition and Information Processing (PRIP) (pp. 140–145).
% -------------------------------------------------------------------------
% INPUTS:
% - ROIonly: Smallest box containing the ROI, with the imaging data ready 
%            for texture analysis computations. Voxels outside the ROI are 
%            set to NaNs.
% - levels: Vector containing the quantized gray-levels in the tumor region
%           (or reconstruction levels of quantization).
%
% ** 'ROIonly' and 'levels' should be outputs from 'prepareVolume.m' **
% -------------------------------------------------------------------------
% OUTPUTS:
% - GLSZM: Gray-Level Size Zone Matrix of 'ROIOnly'.
% -------------------------------------------------------------------------
% AUTHOR(S): Martin Vallieres <mart.vallieres@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: January 2013
% - Revision: May 2015
% -------------------------------------------------------------------------
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


% PRELIMINARY
nLevel = length(levels);
if nLevel > 100
    adjust = 10000;
else
    adjust = 1000;
end
levelTemp = max(levels) + 1;
ROIOnly(isnan(ROIOnly)) = levelTemp;
levels = [levels,levelTemp];


% QUANTIZATION EFFECTS CORRECTION
% In case (for example) we initially wanted to have 64 levels, but due to
% quantization, only 60 resulted.
uniqueVect = round(levels*adjust)/adjust;
ROIOnly = round(ROIOnly*adjust)/adjust;
NL = length(levels) - 1;


% INITIALIZATION
nInit = numel(ROIOnly);
GLSZM = zeros(NL,nInit);


% COMPUTATION OF GLSZM
temp = ROIOnly;
for i = 1:NL
    temp(ROIOnly~=uniqueVect(i)) = 0;
    temp(ROIOnly==uniqueVect(i)) = 1;
    connObjects = bwconncomp(temp,26);
    nZone = length(connObjects.PixelIdxList);
    for j = 1:nZone
        col = length(connObjects.PixelIdxList{j});
        GLSZM(i,col) = GLSZM(i,col) + 1;
    end
end


% REMOVE UNECESSARY COLUMNS
stop = find(sum(GLSZM),1,'last');
GLSZM(:,(stop+1):end) = [];

end