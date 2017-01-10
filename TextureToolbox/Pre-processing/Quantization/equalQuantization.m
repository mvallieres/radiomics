function [ROIonlyEqual,levels] = equalQuantization(ROIonly,Ng)
% -------------------------------------------------------------------------
% function [ROIonlyEqual,levels] = equalQuantization(ROIonly,Ng)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function computes Equal-probability quantization on the region of 
% interest (ROI) of an input volume..
% -------------------------------------------------------------------------
% INPUTS:
% - ROIonly: 3D array, with voxels outside the ROI set to NaNs.
% - Ng: Integer specifying the number of gray levels in the quantization.
% -------------------------------------------------------------------------
% OUTPUTS:
% - ROIonlyEqual: Quantized input volume.
% - levels: Vector containing the quantized gray-level values in the ROI
%           (or reconstruction levels of quantization).
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

ROIonly = double(ROIonly);
ROIonly = (ROIonly - min(ROIonly(:)))./(max(ROIonly(:)) - min(ROIonly(:)));
qmax = max(ROIonly(:));
qmin = min(ROIonly(:));
ROIonly(isnan(ROIonly)) = qmax + 1;

[b,perm] = sort(reshape(ROIonly,[1 size(ROIonly,1)*size(ROIonly,2)*size(ROIonly,3)]));

ind = find(b<=1);
d = b(ind);

% EQUALIZATION USING MATLAB. THIS FUNCTION FORCES THE TRANSFORMATION
% TO BE MONOTONIC. THIS PROPERTY IS DESIRABLE TO PRESERVE INTRINSIC TEXTURES.
J = histeq(d,Ng);

b(ind(1:end)) = J(1:end);
ROIonlyEqual = zeros(1,size(ROIonly,1)*size(ROIonly,2)*size(ROIonly,3));
ROIonlyEqual(perm(1:end)) = b(1:end);

ROIonlyEqual = reshape(ROIonlyEqual,[size(ROIonly,1),size(ROIonly,2),size(ROIonly,3)]);
ROIonlyEqual(ROIonly==2) = NaN;

levels = unique(J);
volumeTemp = ROIonlyEqual;
for i=1:numel(levels)
    ROIonlyEqual(volumeTemp==levels(i)) = i;
end

levels = 1:Ng;

end