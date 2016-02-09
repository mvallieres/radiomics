function [textures] = getGlobalTextures(ROIonly,Nbins)
% -------------------------------------------------------------------------
% function [textures] = getGlobalTextures(ROIonly,Nbins)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function computes Global texture features from the region of 
% interest (ROI) of an input volume.
% -------------------------------------------------------------------------
% INPUTS:
% - ROIonly: Smallest box containing the ROI, with the imaging data ready 
%            for texture analysis computations. Voxels outside the ROI are 
%            set to NaNs.
% - Nbins: Numerical value specifying the number of bins to use in the
%          histogram prior to global textures calculations.
%
% ** 'ROIonly' should be output from 'prepareVolume.m' **
% -------------------------------------------------------------------------
% OUTPUTS:
% - textures: Struture specifying the values of different Global texture
%             features as defined below.
% -------------------------------------------------------------------------
% AUTHOR(S): Martin Vallieres <mart.vallieres@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: January 2013.
% - Revision: May 2015
% -------------------------------------------------------------------------
% STATEMENT:
% This file is part of <https://github.com/mvallieres/radiomics/>, 
% a package providing MATLAB programming tools for radiomics analysis.
% --> Copyright (C) 2013-2015  Martin Vallieres
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


% PRELIMINARY
vectorValid = ROIonly(~isnan(ROIonly));
histo = hist(vectorValid,Nbins);
histo = histo./(sum(histo(:)));
vectNg = 1:Nbins;
u = histo*vectNg';


% COMPUTATION OF TEXTURES
% 1. Variance
variance = 0;
for i=1:Nbins
    variance = variance+histo(i)*(i-u)^2;
end
sigma = sqrt(variance);
textures.Variance = variance;

% 2. Skewness
skewness = 0;
for i = 1:Nbins
    skewness = skewness+histo(i)*(i-u)^3;
end
skewness = skewness/sigma^3;
textures.Skewness = skewness;

% 3. Kurtosis
kurtosis = 0;
for i = 1:Nbins
    kurtosis = kurtosis+histo(i)*(i-u)^4;
end
kurtosis = (kurtosis/sigma^4) - 3;
textures.Kurtosis = kurtosis;

end
