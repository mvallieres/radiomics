function [sDataCT] = correctPatient30_HN(pathDATA)
% -------------------------------------------------------------------------
% function [sDataCT] = correctPatient30_HN(pathDATA)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% Reason of existence of the function: the RTstruct of the CT scan of
% patient 30 is incomplete (most probably it was previously incorrectly 
% processed by MIM). The problem is solved by copying and processing the 
% 'contour' structure of the PET scan.
% -------------------------------------------------------------------------
% INPUTS:
% - pathDATA: Full path to the 'DATA' folder of the HN workspace.
% -------------------------------------------------------------------------
% OUTPUTS:
% - sData: Corrected sData.
% -------------------------------------------------------------------------
% AUTHOR(S): Martin Vallieres <mart.vallieres@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: July 2015
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
cd(pathDATA)
load('Patient30_CT'), sDataCT = sData;
load('Patient30_PET'), sDataPET = sData;

% ADJUSTING CT CONTOURS ACCORDING TO PET CONTOURS.
sDataCT{2}.scan.contour = sDataPET{2}.scan.contour;
nContour = length(sDataPET{2}.scan.contour);
downF = sDataCT{2}.scan.pixelW/sDataPET{2}.scan.pixelW;
offsetStart = round(abs(sDataPET{3}(1).ImagePositionPatient(2)-sDataCT{3}(1).ImagePositionPatient(2))*sDataCT{2}.scan.pixelW);
offsetEnd = round((size(sDataCT{2}.scan.volume,1)*sDataCT{2}.scan.pixelW - size(sDataPET{2}.scan.volume,1)*sDataPET{2}.scan.pixelW)/sDataCT{2}.scan.pixelW)-offsetStart;
for i = 1:nContour
    sDataCT{2}.scan.contour(i).boxBound(1:2,1) = round((sDataPET{2}.scan.contour(i).boxBound(1:2,1)/downF + offsetStart)); sDataCT{2}.scan.contour(i).boxBound(1:2,2) = round((sDataPET{2}.scan.contour(i).boxBound(1:2,2)/downF + offsetEnd));
    boxBound = sDataCT{2}.scan.contour(i).boxBound;
    nSlices = size(sDataCT{2}.scan.contour(i).boxMask,3);
    sDataCT{2}.scan.contour(i).boxMask = zeros(boxBound(1,2)-boxBound(1,1)+1,boxBound(2,2)-boxBound(2,1)+1,nSlices);
    for j = 1:nSlices
        sDataCT{2}.scan.contour(i).boxMask(:,:,j) = imresize(sDataPET{2}.scan.contour(i).boxMask(:,:,j),[boxBound(1,2)-boxBound(1,1)+1,boxBound(2,2)-boxBound(2,1)+1],'nearest');
    end
end

cd(startpath)
end