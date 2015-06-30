function [volumeWBPF] = waveletBPfilt(volume,ratio,wavelet)
% -------------------------------------------------------------------------
% [volumeWBPF] = waveletBPfilt(volume,ratio,wavelet)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function performs wavelet band-pass filtering to the input volume.
% -------------------------------------------------------------------------
% REFERENCE:
% [1] Vallieres, M. et al. (2015). A radiomics model from joint FDG-PET and 
%     MRI texture features for the prediction of lung metastases in soft-tissue 
%     sarcomas of the extremities. Physics in Medicine and Biology, 60(14), 
%     5471-5496. doi:10.1088/0031-9155/60/14/5471
% -------------------------------------------------------------------------
% INPUTS:
% - volume: 3D array containing the volume to filter
% - ratio: Numerical value specifying the 'R' ratio as defined in ref. [1].
% - wavelet: String specifying the name of the MATLAB wavelet to use. 
%            Example: 'sym8'
% -------------------------------------------------------------------------
% OUTPUTS:
% - volumeWBPF: Filtered input volume.
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

weightHHH_LLL=8/(2*(3*ratio+1));
weight_Rest=8*ratio/(2*(3*ratio+1));

wDEC=wavedec3(volume,1,wavelet);
nbcell=length(wDEC.dec);

wDEC.dec{1}=weightHHH_LLL.*wDEC.dec{1};
for i=2:nbcell-1
   wDEC.dec{i}=weight_Rest.*wDEC.dec{i};
end
wDEC.dec{end}=weightHHH_LLL.*wDEC.dec{end};

volumeWBPF=waverec3(wDEC);

end