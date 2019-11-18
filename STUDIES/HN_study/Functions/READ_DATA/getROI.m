function [ROIbox,ROImask,boxBound] = getROI(sData,contourNumber)
% -------------------------------------------------------------------------
% function [ROIbox,ROImask,boxBound] = getROI(sData,contourNumber)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% Computes the ROI box (smallest box containing the region of interest) 
% and associated mask from a 'sData' file.
% -------------------------------------------------------------------------
% INPUTS:
% - sData: Cell of structures organizing the data.
% - contourNumber: Which contour to use in the computation of the tumor
%                  box. See sData{2}.scan.contour(contourNumber).
% -------------------------------------------------------------------------
% OUTPUTS:
% - ROIbox: 3D array of imaging data defining the smallest box
%           containing the region of interest.
% - ROImask: 3D array of 1's and 0's defining the ROI in ROIbox.
% - boxBound: Bounds of the smallest box containing the ROI in the volume.
%             Format: [minRow, maxRow;
%                      minColumn, maxColumns;
%                      minSlice, maxSlice]
% -------------------------------------------------------------------------
% AUTHOR(S): 
% - Martin Vallieres <mart.vallieres@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: February 2016
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

% GETTING THE XYZ POINTS FROM THE sData STRUCTURE
ROI_XYZ = sData{2}.scan.contour(contourNumber).points_XYZ;
spatialRef = sData{2}.scan.volume.spatialRef;
ROIbox = sData{2}.scan.volume.data;


% COMPUTING THE ROI
sz = size(ROIbox);
ROImask = zeros(sz(1),sz(2),sz(3));
[X,Y,Z] = worldToIntrinsic(spatialRef,ROI_XYZ(:,1),ROI_XYZ(:,2),ROI_XYZ(:,3)); % X,Y,Z in intrinsic image coordinates
points = [X,Y,Z];
if strcmp(sData{2}.scan.orientation,'Axial')
    a = 1; b = 2; c = 3;
elseif strcmp(sData{2}.scan.orientation,'Sagittal')
    a = 2; b = 3; c = 1;
elseif strcmp(sData{2}.scan.orientation,'Coronal')
    a = 1; b = 3; c = 2;
end
K = round(points(:,c)); % Must assign the points to one slice
slices = unique(K);
for k = 1:numel(slices)
    ind = find(K == slices(k));
    ROImask(:,:,slices(k)) = or(ROImask(:,:,slices(k)),poly2mask(points(ind,a),points(ind,b),sz(1),sz(2)));
end


% COMPUTING THE SMALLEST BOX CONTAINING THE REGION OF INTEREST
[boxBound] = computeBoundingBox(ROImask);
ROIbox = ROIbox(boxBound(1,1):boxBound(1,2),boxBound(2,1):boxBound(2,2),boxBound(3,1):boxBound(3,2));
ROImask = ROImask(boxBound(1,1):boxBound(1,2),boxBound(2,1):boxBound(2,2),boxBound(3,1):boxBound(3,2));

end
