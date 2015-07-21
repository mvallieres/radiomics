function [percentInactive] = getPercentInactive(ROIonlyPET,thresh)
% -------------------------------------------------------------------------
% function [percentInactive] = getPercentInactive(ROIonlyPET,thresh)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function computes the percentage of inactive PET volume from the 
% region of interest (ROI) of an input volume. A typical threshold of 
% thresh × (SUVmax)^2 followed by closing and opening morphological 
% operations is used to differentiate active and inactive regions.
% -------------------------------------------------------------------------
% INPUTS:
% - ROIonlyPET: 3D array representing the PET volume in SUV format, with 
%               voxels outside the ROI set to NaNs.
% - thresh: Numerical value specifying the threshold threshold thresh × 
%           (SUVmax)^2.
% -------------------------------------------------------------------------
% OUTPUTS:
% - percentInactive: Percentage of the ROI that is inactive.
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

% APPLYING INITIAL THRESHOLD
mask = ~isnan(ROIonlyPET);
ROIonlyPET(isnan(ROIonlyPET)) = 0;
maskInactive = ROIonlyPET > thresh*(max(ROIonlyPET(:)))^2;


% MORPHOLOGICAL OPERATIONS
conn = zeros(5,5,5);
conn(:,:,1) = [0 0 0 0 0;0 0 0 0 0;0 0 1 0 0;0 0 0 0 0;0 0 0 0 0];
conn(:,:,5) = conn(:,:,1);
conn(:,:,2) = [0 0 0 0 0;0 0 1 0 0;0 1 1 1 0;0 0 1 0 0;0 0 0 0 0];
conn(:,:,4) = conn(:,:,2);
conn(:,:,3) = [0 0 1 0 0;0 1 1 1 0;1 1 1 1 1;0 1 1 1 0;0 0 1 0 0];
maskInactive = imclose(maskInactive,conn);
maskInactive = imopen(maskInactive,conn);


% COMPUTING THE PERCENTAGE OF INACTIVE VOLUME
perimeter = bwperim(mask,26);
newMask = maskInactive+perimeter;
newMask(mask==0)=10; % temporary value
newMask(newMask==1)=10;
newMask(newMask==2)=10;
newMask(newMask==0) = 1;
newMask(newMask==10) = 0;
connObjects = bwconncomp(newMask,26);
b = zeros(1,connObjects.NumObjects);
for i = 1:connObjects.NumObjects
    a = find(length(connObjects.PixelIdxList{i})>=15); % If the number of of pixel forming and object is lower than 15, reject it
    if isempty(a)
        b(i) = 0;
    else
        b(i) = 1;
    end
end
[row,col] = find(b>0);
sumInactive = 0;
for i = 1:length(col)
    sumInactive = sumInactive + length(connObjects.PixelIdxList{row(i),col(i)});
end
sumVolume = sum(mask(:));
percentInactive = sumInactive/sumVolume*100;

end