function [ROIonlyLloyd,levels] = lloydQuantization(ROIonly,Ng)
% -------------------------------------------------------------------------
% function [ROIonlyLloyd,levels] = lloydQuantization(ROIonly,Ng)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function computes Lloyd-Max quantization on the region of interest 
% (ROI) of an input volume.
% -------------------------------------------------------------------------
% REFERENCES:
% [1] Max, J. (1960). Quantizing for minimum distortion. IEEE Trans. Inf.
%     Theory, 6(1), 7–12.
% [2] Lloyd, S. (1982). Least Squanres Quantization in PCM. IEEE Trans.
%     Inf. Theory, IT-28(2),129–137.
% -------------------------------------------------------------------------
% INPUTS:
% - ROIonly: 3D array, with voxels outside the ROI set to NaNs.
% - Ng: Integer specifying the number of gray levels in the quantization.
% -------------------------------------------------------------------------
% OUTPUTS:
% - ROIonlyLloyd: Quantized input volume.
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

d = b(find(b<=1));

[transitions,~] = lloyds(d,Ng); % LLoyd- algorithm

temp = [qmin,transitions,qmax];
transitions = temp;

sortEqual = zeros(1,size(ROIonly,1)*size(ROIonly,2)*size(ROIonly,3));

sortEqual(find(b <= transitions(2) & b >= transitions(1))) = 1;
for i = 2:Ng
    sortEqual(find(b <= transitions(i+1) & b > transitions(i))) = i;
end

ROIonlyLloyd = zeros(1,size(ROIonly,1)*size(ROIonly,2)*size(ROIonly,3));
ROIonlyLloyd(perm(1:end)) = sortEqual(1:end);

ROIonlyLloyd = reshape(ROIonlyLloyd,[size(ROIonly,1),size(ROIonly,2),size(ROIonly,3)]);
ROIonlyLloyd(ROIonly==2) = NaN;

levels = 1:Ng;

end