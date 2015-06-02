function [boxBound] = computeBoundingBox(mask)
% -------------------------------------------------------------------------
% function [boxBound] = computeBoundingBox(mask)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function computes the smallest box containing the whole region of 
% interest (ROI). It is adapted from the function compute_boundingbox.m
% of CERR <http://www.cerr.info/>.
% -------------------------------------------------------------------------
% INPUTS:
% - mask: 3D array, with 1's inside the ROI, and 0's outside the ROI.
% -------------------------------------------------------------------------
% OUTPUTS:
% - boxBound: Bounds of the smallest box containing the ROI. 
%             Format: [minRow, maxRow;
%                      minColumn, maxColumns;
%                      minSlice, maxSlice]
% -------------------------------------------------------------------------
% AUTHOR(S): 
% - Martin Vallieres <mart.vallieres@gmail.com>
% - CERR development team <http://www.cerr.info/>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: May 2015
%--------------------------------------------------------------------------
% STATEMENT:
% This file is part of <https://github.com/mvallieres/radiomics/>, 
% a package providing MATLAB programming tools for radiomics analysis.
% --> Copyright (C) 2015  Martin Vallieres
% --> Copyright 2010, Joseph O. Deasy, on behalf of the CERR development team
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

[iV,jV,kV] = find3d(mask);
boxBound(1,1) = min(iV);
boxBound(1,2) = max(iV);
boxBound(2,1) = min(jV);
boxBound(2,2) = max(jV);
boxBound(3,1) = min(kV);
boxBound(3,2) = max(kV);

end


% CERR UTILITY FUNCTIONS (can be found at: https://github.com/adityaapte/CERR)
function [iV,jV,kV] = find3d(mask3M)
indV = find(mask3M(:));
[iV,jV,kV] = fastind2sub(size(mask3M),indV);
iV = iV';
jV = jV';
kV = kV';
end

function varargout = fastind2sub(siz,ndx)
nout = max(nargout,1);
if length(siz)<=nout,
  siz = [siz ones(1,nout-length(siz))];
else
  siz = [siz(1:nout-1) prod(siz(nout:end))];
end
n = length(siz);
k = [1 cumprod(siz(1:end-1))];
ndx = ndx - 1;
for i = n:-1:1,
  varargout{i} = floor(ndx/k(i)) + 1;
  ndx = ndx - (varargout{i}-1) * k(i);
end
end