function [NGTDM,countValid] = getNGTDM(ROIOnly,levels)
% -------------------------------------------------------------------------
% function [NGTDM,countValid] = getNGTDM(ROIOnly,levels)
% -------------------------------------------------------------------------
% DESCRIPTION:
% This function computes the Neighborhood Gray-Tone Difference Matrix 
% (NGTDM) of the region of interest (ROI) of an input volume. The input 
% volume is assumed to be isotropically resampled. The NGTDM is computed 
% using 26-voxel connectivity. To account for discretization length 
% differences, all averages around a center voxel are performed such that 
% the neighbours at a distance of sqrt(3) voxels are given a weight of 
% 1/sqrt(3), and the neighbours at a distance of sqrt(2) voxels are given a 
% weight of 1/sqrt(2).
%
% --> This function is compatible with 2D analysis (language not adapted in the text)
% -------------------------------------------------------------------------
% REFERENCE:
% [1] Amadasun, M., & King, R. (1989). Textural Features Corresponding to 
%     Textural Properties. IEEE Transactions on Systems Man and Cybernetics,
%     19(5), 1264–1274.
% -------------------------------------------------------------------------
% INPUTS:
% - ROIonly: Smallest box containing the ROI, with the imaging data ready 
%            for texture analysis computations. Voxels outside the ROI are 
%            set to NaNs.
% - levels: Vector containing the quantized gray-levels in the tumor region
%           (or reconstruction levels of quantization).
%
% ** 'ROIonly' and 'levels' should be outputs from 'prepareVolume.m' **
% -------------------------------------------------------------------------
% OUTPUTS:
% - NGTDM: Neighborhood Gray-Tone Difference Matrix of 'ROIOnly'.
% - countValid: Number of valid voxels used in the NGTDM computation. 
%               Required for the computation of texture features in 
%               'getNGTDMtextures.m'
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


% PRELIMINARY
if numel(size(ROIOnly)) == 2 % generalization to 2D inputs
    twoD = 1;
else
    twoD = 0;
end
nLevel = length(levels);
if nLevel > 100
    adjust = 10000;
else
    adjust = 1000;
end
if twoD
    ROIOnly = padarray(ROIOnly,[1,1],NaN,'both');
else
    ROIOnly = padarray(ROIOnly,[1,1,1],NaN,'both');
end


% QUANTIZATION EFFECTS CORRECTION
% In case (for example) we initially wanted to have 64 levels, but due to
% quantization, only 60 resulted.
uniqueVol = round(levels*adjust)/adjust;
ROIOnly = round(ROIOnly*adjust)/adjust;
NL = length(levels);
temp = ROIOnly;
for i = 1:NL
    ROIOnly(temp==uniqueVol(i)) = i;
end


% INTIALIZATION
NGTDM = zeros(NL,1);
countValid = zeros(NL,1);


% COMPUTATION OF NGTDM
if twoD
    [i,j] = ind2sub(size(ROIOnly),find(~isnan(ROIOnly)));
    posValid = [i,j];
    nValid_temp = size(posValid,1);
    w8 = 1/sqrt(2); w4 = 1; % Weights given to different neighbors to correct for discretization length differences
    weights = [w8 w4 w8; w4 w4 w4; w8 w4 w8];
    weights = weights(:);
    for n = 1:nValid_temp
        neighbours = zeros(9,1);
        neighbours(1:9) = ROIOnly((posValid(n,1)-1):(posValid(n,1)+1),(posValid(n,2)-1):(posValid(n,2)+1));
        neighbours = neighbours.*weights;
        value = neighbours(5);
        neighbours(5) = NaN;
        neighbours = neighbours/sum(weights(~isnan(neighbours)));
        neighbours(5) = []; % Remove the center voxel
        if ~isempty(neighbours(~isnan(neighbours))) % Thus only excluding voxels with NaNs only as neighbors.
            NGTDM(value) = NGTDM(value) + abs(value-sum(neighbours(~isnan(neighbours))));
            countValid(value) = countValid(value) + 1;
        end
    end
else
    [i,j,k] = ind2sub(size(ROIOnly),find(~isnan(ROIOnly)));
    posValid = [i,j,k];
    nValid_temp = size(posValid,1);
    w26 = 1/sqrt(3); w18 = 1/sqrt(2); w6=1; % Weights given to different neighbors to correct for discretization length differences
    weights = [w26 w18 w26; w18 w6 w18; w26 w18 w26];
    weights(:,:,2) = [w18 w6 w18; w6 w6 w6; w18 w6 w18];
    weights(:,:,3) = weights(:,:,1);
    weights = weights(:);
    for n = 1:nValid_temp
        neighbours = zeros(27,1);
        neighbours(1:27) = ROIOnly((posValid(n,1)-1):(posValid(n,1)+1),(posValid(n,2)-1):(posValid(n,2)+1),(posValid(n,3)-1):(posValid(n,3)+1));
        neighbours = neighbours.*weights;
        value = neighbours(14);
        neighbours(14) = NaN;
        neighbours = neighbours/sum(weights(~isnan(neighbours)));
        neighbours(14) = []; % Remove the center voxel
        if ~isempty(neighbours(~isnan(neighbours))) % Thus only excluding voxels with NaNs only as neighbors.
            NGTDM(value) = NGTDM(value) + abs(value-sum(neighbours(~isnan(neighbours))));
            countValid(value) = countValid(value) + 1;
        end
    end
end

end