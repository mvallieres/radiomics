function [fusedBox] = fusePETCTissam(ROIbox_PET,ROIbox_CT,maskBox_PET,CTinv,CTweight,wavelet,fusWeights)
% -------------------------------------------------------------------------
% function [fusedBox] = fusePETCT(ROIbox_PET,ROIbox_CT,maskBox_PET,CTinv,CTweight,wavelet)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function fuses the region of interest (ROI) of two registered PET and 
% CT volumes using a technique based on the wavelet transform. See Ref. [2] 
% for more details.
% -------------------------------------------------------------------------
% REFERENCE:
% [1] Vallieres, M. et al. (2015). FDG-PET/CT radiomics models for the 
%     early prediction of different tumour outcomes in head and neck cancer.
%     The Journal of Nuclear Medicine, aa(bb), xxx-yyy. 
%     doi:
% -------------------------------------------------------------------------
% INPUTS:
% - ROIbox_PET: 3D array of the smallest box containing the ROI of the PET 
%               volume.
% - ROIbox_CT: 3D array of the smallest box containing the ROI of the CT
%               volume.
% - maskBox_PET: 3D array of the mask of the smallest box specifying the
%                ROI of the PET volume. Voxels within the ROI are assigned
%                value of 1, voxels outside a value of 0.
% - CTinv: String specifying if the intensities of the CT volume are
%          inverted prior to fusion with PET. Either 'Inv' for inversion,
%          or 'NoInv' for no inversion.
% - CTweight: Numerical value specifying the weight of the CT scan in the
%              fusion with PET.
% - wavelet: String specifying the name of the MATLAB wavelet basis used
%            wavelet basis used in the fusion process.
% -------------------------------------------------------------------------
% OUTPUTS:
% - fusedBox: 3D array of the smallest box containing the ROI of the fused
%             PET/CT volume.
% -------------------------------------------------------------------------
% AUTHOR(S): Martin Vallieres <mart.vallieres@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: July 2015
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


% PET PRE_PROCESSING
ROIbox_PET = sqrt(ROIbox_PET);


% RESAMPLING CT VOLUME TO PET IN-PLANE RESOLUTION
szPET = size(ROIbox_PET);
temp=zeros(szPET);
for k = 1:szPET(3)
    temp(:,:,k) = imresize(ROIbox_CT(:,:,k),[szPET(1),szPET(2)],'Method','cubic','Antialiasing',true);
end
ROIbox_CT = temp;


% NORMALIZATION (and inversion) OF VOLUMES
ROIOnly_CT = ROIbox_CT;
ROIOnly_PET = ROIbox_PET;
ROIOnly_CT(maskBox_PET==0) = NaN;
ROIOnly_PET(maskBox_PET==0) = NaN;
minCT = min(ROIOnly_CT(:));
ROIbox_CT = ROIbox_CT - minCT;
ROIOnly_CT = ROIOnly_CT - minCT;
ROIbox_CT = ROIbox_CT./max(ROIOnly_CT(:)).*255;
if strcmp(CTinv,'Inv')
    ROIbox_CT = 255 - ROIbox_CT;
end
minPET = min(ROIOnly_PET(:));
ROIbox_PET = ROIbox_PET - minPET;
ROIOnly_PET = ROIOnly_PET - minPET;
ROIbox_PET = ROIbox_PET./max(ROIOnly_PET(:)).*255;


% WAVELET DECOMPOSITION AND FUSION OF COEFFICIENTS
wdecCT = wavedec3(ROIbox_CT,1,wavelet);
wdecPET = wavedec3(ROIbox_PET,1,wavelet);
wdecFUSED = wdecPET;
szDEC = size(wdecFUSED.dec{1}); % All sub-bands are of the same size

% Fusion of the LLL sub-bands (dec{1})
b = 1;
wdecFUSED.dec{b} = (1-CTweight).*wdecPET.dec{b} + CTweight.*wdecCT.dec{b};

% Fusion of the rest of ths sub-bands
for b = 2:8
    wdecFUSED.dec{b} = fusWeights((b-2)*2 + 1).*wdecPET.dec{b} + fusWeights((b-1)*2).*wdecCT.dec{b};
end

% WAVELET RECONSTRUCTION
fusedBox = waverec3(wdecFUSED);
fusedOnly = fusedBox;
fusedOnly(maskBox_PET==0) = NaN;
minFus = min(fusedOnly(:));
fusedBox = fusedBox - minFus;
fusedOnly = fusedOnly - minFus;
fusedBox = fusedBox./max(fusedOnly(:)).*255;

end