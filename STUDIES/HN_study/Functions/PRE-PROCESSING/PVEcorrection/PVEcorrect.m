function [volumePVC] = PVEcorrect(volume,nIter,waveletName)
% -------------------------------------------------------------------------
% function [volumePVC] = PVEcorrect(volume,nIter,waveletName)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% Apply partial-volume effect (PVE) correction of an input PET volume
% using the methodology developed in ref. [1].
% -------------------------------------------------------------------------
% REFERENCE:
% [1] Boussion, N. et al. (2009). Incorporation of wavelet-based denoising
%     in iterative deconvolution for partial volume correction in 
%     whole-body PET imaging. Eur J Nucl Med Mol Imaging, 36(7), 1064-1075.
% -------------------------------------------------------------------------
% INPUTS:
% - volume: 3D array representing the input PET volume to correct for PVE.
% - nIter: Number of iterations in the deconvolution process of the PVE 
%          correction (see ref. [1]). If set to a string 'compute', the 
%          algorithm will determine the optimal number of iterations in 
%          terms of residuals up to a maximum of 10.
% - waveletName: (optional). MATLAB name of the type of wavelet used in the 
%                denoising part of the PVE correction (see ref. [1]). 
%                Default is 'bior3.5'.
% -------------------------------------------------------------------------
% OUTPUTS:
% - volumePVC: 3D array representing the input PET volume corrected for
%              partial volume effects.
% -------------------------------------------------------------------------
% AUTHOR(S): 
% - Andre Diamant <adboustead@gmail.com>
% - Martin Vallieres <mart.vallieres@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: January 2016
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

% INTIALIZATION
volumePVC = zeros(size(volume));
psf = nonIsotropicGaussianPSF(1,4);

if nargin == 3
    global wavelet_name
    wavelet_name = waveletName;
end

% PVE CORRECTION
if strcmp(nIter,'compute')
    [~,residuals] = deconvlucydenoiseNew(volume,psf,10);
    [~,nIter] = min(abs(residuals));
    fprintf('*** FOUND NUMBER OF ITERATIONS TO BE %.0f ***\n',nIter)
end
[volumePVC,~] = deconvlucydenoiseNew(volume,psf,nIter);

if nargin == 3
    clear global wavelet_name
end

end