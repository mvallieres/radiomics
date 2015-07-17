function [solidity] = getSolidity(ROIonly,pixelW,sliceS)
% -------------------------------------------------------------------------
% function [solidity] = getSolidity(ROIonly,pixelW,sliceS)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function computes the solidity metric of the region of interest 
% (ROI) of an input volume.
% -------------------------------------------------------------------------
% INPUTS:
% - ROIonly: 3D array, with voxels outside the ROI set to NaNs.
% - pixelW: Pixel width, or in-plane resolution, in mm.
% - sliceS: Slice spacing, in mm.
% -------------------------------------------------------------------------
% OUTPUTS:
% - solidity: Ratio of the number of voxels in the ROI to the number of 
%             voxels in the 3D convex hull of the ROI (smallest polyhedron 
%             containing the ROI).
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

mask = ~isnan(ROIonly); % Find mask covering the ROI

% ISOTROPIC RESAMPLING
sFactor = sliceS/pixelW; % scaling factor
mask = imresize3D(mask,[],[round(double(size(mask,1))),round(double(size(mask,2))),round(double(size(mask,3))*sFactor)],'nearest','fill');

% SOLIDITY COMPUTATION
perimeter = bwperim(mask,18);
nPoints = length(find(perimeter));
X = zeros(nPoints,1); Y = zeros(nPoints,1); Z = zeros(nPoints,1);
count = 1;
for i = 1:size(perimeter,3)
    [row,col] = find(perimeter(:,:,i));
    p = length(row);
    if p > 0
        X(count:count+p-1,1) = col(1:end);
        Y(count:count+p-1,1) = row(1:end);
        Z(count:count+p-1,1) = i;
        count = count + p;
    end
end
[~,volumeH] = convhull(X,Y,Z);
volumeROI = sum(mask(:));
solidity = volumeROI/volumeH;

end