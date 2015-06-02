function [ROIonlyNorm] = CollewetNorm(ROIonly)
% -------------------------------------------------------------------------
% function [ROIonlyNorm] = CollewetNorm(ROIonly)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function applies the intensity normalization scheme found in 
% reference [1]. It basically identifies voxels of the region of interest 
% (ROI) of an input volume with intensities outside mean + 3*std, and 
% equals them to NaNs.
% -------------------------------------------------------------------------
% REFERENCE:
% [1] Collewet, G., Strzelecki, M. and Mariette F. (2004). Influence of MRI 
%     acquisition protocols and image intensity normalization methods on 
%     texture classification. Magn. Reson. Imaging, 22(1), 81–91.
% -------------------------------------------------------------------------
% INPUTS:
% - ROIonly: 3D array, with voxels outside the ROI set to NaNs.
% -------------------------------------------------------------------------
% OUTPUTS:
% - ROIonlyNorm: Normalized input volume using Collewet method.
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
temp = ROIonly(~isnan(ROIonly));
u = mean(temp);
sigma = std(temp);
ROIonlyNorm = ROIonly;
ROIonlyNorm(ROIonly > (u + 3*sigma)) = NaN;
ROIonlyNorm(ROIonly < (u - 3*sigma)) = NaN;

end