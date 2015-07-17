function [volume] = getVolume(ROIonly,pixelW,sliceS)
% -------------------------------------------------------------------------
% function [volume] = getVolume(ROIonly,pixelW,sliceS)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function computes the volume an input region region of interest (ROI).
% -------------------------------------------------------------------------
% INPUTS:
% - ROIonly: 3D array, with voxels outside the ROI set to NaNs.
% - pixelW: Pixel width, or in-plane resolution, in mm.
% - sliceS: Slice spacing, in mm.
% -------------------------------------------------------------------------
% OUTPUTS:
% - volume: Number of voxels in the ROI multiplied by the voxel dimension, 
%           in mm^3.
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
numberVoxel = sum(mask(:));
volume = numberVoxel * pixelW * pixelW * sliceS;

end