function [ROIonly] = getROIonly(sData,contourNumber)
% -------------------------------------------------------------------------
% function [ROIonly] = getROIonly(sData,contourNumber)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% Computes the ROI box (smallest box containing the region of interest) 
% from a 'sData' file. Then, voxels outside the tumor region are set to NaNs.
% -------------------------------------------------------------------------
% INPUTS:
% - sData: Cell of structures organizing the data.
% - contourNumber: Which contour to use in the computation of the tumor
%                  box. See sData{2}.scan.contour(contourNumber). 
% -------------------------------------------------------------------------
% OUTPUTS:
% - ROIonly: 3D array of imaging data defining the smallest box containing 
%            the region of interest. The voxels outside the tumor region 
%            are set to NaNs.
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


ROIonly = getROIbox(sData,contourNumber);
mask = sData{2}.scan.contour(contourNumber).boxMask;
ROIonly(~mask) = NaN;

end
