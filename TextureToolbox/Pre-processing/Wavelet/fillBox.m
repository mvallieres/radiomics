function [ROIboxFill] = fillBox(ROIonly)
% -------------------------------------------------------------------------
% function [ROIboxFill] = fillBox(ROIonly)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function randomly fills in the NaNs outside the region of interest 
% (ROI) of the input volume according to its intensity distribution.
% -------------------------------------------------------------------------
% INPUTS:
% - ROIonly: 3D array, with voxels outside the ROI set to NaNs.
% -------------------------------------------------------------------------
% OUTPUTS:
% - ROIboxFill: Input volume for which NaNs have been filled in according 
%               to the intensity distribution of the ROI.
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


% INITIALIZATION
nBins = 128; % Number of bins in the histogram, default value
nFill = sum(sum(sum(isnan(ROIonly)))); % Number of positions to fill
[row,col,z] = ind2sub(size(ROIonly),find(isnan(ROIonly)));
index = [row,col,z];


% BUILDING THE CUMULATIVE FUNCTION (cumFunct)
values = ROIonly(~isnan(ROIonly));
[histo,val] = hist(values,nBins);
histo = histo./(sum(histo(:))); 
cumFunct = cumsum(histo);


% RANDOM NUMBER GENERATOR SEED 
rng(1000000); % Set the random generator to a specific seed, to obtain reproducible fusions.


% FILLING THE BOX
randVect = rand(1,nFill);
ROIboxFill = ROIonly;
for i = 1:nFill
    ROIboxFill(index(i,1),index(i,2),index(i,3)) = val(find(cumFunct > randVect(i),1,'first'));
end

end