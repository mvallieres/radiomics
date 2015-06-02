function [ROIbox] = getROIbox(sData,contourNumber)
% -------------------------------------------------------------------------
% function [ROIbox] = getROIbox(sData,contourNumber)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% Computes the ROI box (smallest box containing the region of interest) 
% from a 'sData' file.
% -------------------------------------------------------------------------
% INPUTS:
% - sData: Cell of structures organizing the data.
% - contourNumber: Which contour to use in the computation of the tumor
%                  box. See sData{2}.scan.contour(contourNumber).
% -------------------------------------------------------------------------
% OUTPUTS:
% - ROIbox: 3D array of imaging data defining the smallest box
%           containing the region of interest.
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


volume = sData{2}.scan.volume;
boxBound = sData{2}.scan.contour(contourNumber).boxBound;
ROIbox = volume(boxBound(1,1):boxBound(1,2), ...
                  boxBound(2,1):boxBound(2,2), ...
                  boxBound(3,1):boxBound(3,2));
              
end
