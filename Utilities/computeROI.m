function [sData] = computeROI(sData)
% -------------------------------------------------------------------------
% function [sData] = computeROI(sData)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function computes the region of interest (ROI) defined by a given 
% RTstruct in a 'sData' file.
% -------------------------------------------------------------------------
% INPUTS: 'sData' file
% -------------------------------------------------------------------------
% OUTPUTS: 'sData' file with ROI definition
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


% INITIALIZATION
sizeScan = size(sData{2}.scan.volume);
sData{2}.scan.contour.name     = [];
sData{2}.scan.contour.boxBound = [];
sData{2}.scan.contour.boxMask  = [];
if strcmp(sData{2}.scan.orientation,'Sagittal') 
    a = 2;
    b = 3;
elseif strcmp(sData{2}.scan.orientation,'Coronal')
    a = 1;
    b = 3;
else 
    a = 1;
    b = 2;
end


% ROI COMPUTATION
nContours = length(fieldnames(sData{4}.ROIContourSequence));
for contourNumber = 1:nContours
    
    % Initialization
    mask = zeros(sizeScan);
    itemContour = ['Item_',num2str(contourNumber)];
    sData{2}.scan.contour(contourNumber).name = sData{4}.StructureSetROISequence.(itemContour).ROIName;
    nSlices = length(fieldnames(sData{4}.ROIContourSequence.(itemContour).ContourSequence));
    
    for sliceNumber = 1:nSlices
        
        % Find slice correspondence between volume and RTstruct
        itemSlice = ['Item_',num2str(sliceNumber)];
        UIDrt = sData{4}.ROIContourSequence.(itemContour).ContourSequence.(itemSlice).ContourImageSequence.Item_1.ReferencedSOPInstanceUID;
        for i = 1:sizeScan(3)
            UIDslice = sData{3}(i).SOPInstanceUID;
            if strcmp(UIDrt,UIDslice)
                sliceOK = i;
                break
            end
        end
        
        pts_temp = sData{4}.ROIContourSequence.(itemContour).ContourSequence.(itemSlice).ContourData; % points stored in the RTstruct file
        if ~isempty(pts_temp)
            
            % Get XYZ points in the reference frame coordinates
            ind = 1:numel(pts_temp)/3;
            pts = zeros([numel(pts_temp)/3,3]);
            pts(:,1) = pts_temp(ind*3-2); pts(:,2) = pts_temp(ind*3-1); pts(:,3) = pts_temp(ind*3);
            pts(:,1) = pts(:,1) - sData{3}(sliceOK).ImagePositionPatient(1);
            pts(:,2) = pts(:,2) - sData{3}(sliceOK).ImagePositionPatient(2);
            pts(:,3) = pts(:,3) - sData{3}(sliceOK).ImagePositionPatient(3);

            % Get transformation matrix
            p1 = sData{3}(sliceOK).PixelSpacing(1);
            p2 = sData{3}(sliceOK).PixelSpacing(2);
            m = [sData{3}(sliceOK).ImageOrientationPatient(a)*p1 sData{3}(sliceOK).ImageOrientationPatient(a+3)*p2; ...
            sData{3}(sliceOK).ImageOrientationPatient(b)*p1 sData{3}(sliceOK).ImageOrientationPatient(b+3)*p2];

            % Transform points from reference frame to image coordinates
            pts = ((m^-1)*pts(:,[a,b])')' + 1; % +1 for MATLAB image coordinates

            % Obtain mask using set of image points
            mask(:,:,sliceOK) = or(mask(:,:,sliceOK),poly2mask(pts(:,1),pts(:,2),sizeScan(1),sizeScan(2)));
        else
            mask(:,:,sliceOK) = zeros(sizeScan(1),sizeScan(2));
        end
    end
    
    % Compute the smallest box containing the whole tumor (and its associated bounds)
    [boxBound] = computeBoundingBox(mask);
    mask = mask(boxBound(1,1):boxBound(1,2),boxBound(2,1):boxBound(2,2),boxBound(3,1):boxBound(3,2));
    sData{2}.scan.contour(contourNumber).boxBound = boxBound;
    sData{2}.scan.contour(contourNumber).boxMask = mask;
end

end