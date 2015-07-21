function [fusedBox] = fusePETMRI(ROIbox_PET,ROIbox_MRI,maskBox_PET,MRIinv,MRIweight,wavelet)
% -------------------------------------------------------------------------
% function [fusedBox] = fusePETMRI(ROIbox_PET,ROIbox_MRI,maskBox_PET,Invert,MRIweight,wavelet)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function fuses the region of interest (ROI) of two registered PET and 
% MRI volumes using a technique based on the wavelet transform. See Ref. [1] 
% for more details.
% -------------------------------------------------------------------------
% REFERENCE:
% [1] Vallieres, M. et al. (2015). A radiomics model from joint FDG-PET and 
%     MRI texture features for the prediction of lung metastases in soft-tissue 
%     sarcomas of the extremities. Physics in Medicine and Biology, 60(14), 
%     5471-5496. doi:10.1088/0031-9155/60/14/5471
% -------------------------------------------------------------------------
% INPUTS:
% - ROIbox_PET: 3D array of the smallest box containing the ROI of the PET 
%               volume.
% - ROIbox_MRI: 3D array of the smallest box containing the ROI of the MRI
%               volume.
% - maskBox_PET: 3D array of the mask of the smallest box specifying the
%                ROI of the PET volume. Voxels within the ROI are assigned
%                value of 1, voxels outside a value of 0.
% - MRIinv: String specifying if the intensities of the MRI volume are
%           inverted prior to fusion with PET. Either 'Inv' for inversion,
%           or 'NoInv' for no inversion.
% - MRIweight: Numerical value specifying the weight of the MRI scan in the
%              fusion with PET.
% - wavelet: String specifying the name of the MATLAB wavelet basis used
%            wavelet basis used in the fusion process.
% -------------------------------------------------------------------------
% OUTPUTS:
% - fusedBox: 3D array of the smallest box containing the ROI of the fused
%             PET/MRI volume. Note that the output contains random values
%             outside the ROI. These values are added in the fusion process 
%             in order to accurately calculate the mean and the standard 
%             deviation of the ROI (only) of the MRI volume when performing
%             Collewet normalization in the wavelet domain. Apply the mask 
%             to 'fusedBox' by setting NaNs outside the ROI to accurately
%             visualize the fusion.
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


% PET PRE_PROCESSING
ROIbox_PET = sqrt(ROIbox_PET);


% RESAMPLING MRI VOLUME TO PET IN-PLANE RESOLUTION (slice spacing
% previously verified to be the same)
szPET = size(ROIbox_PET);
temp=zeros(szPET);
for i = 1:szPET(3)
    temp(:,:,i) = imresize(ROIbox_MRI(:,:,i),[szPET(1),szPET(2)],'Method','cubic','Antialiasing',true);
end
ROIbox_MRI = temp;


% NORMALIZATION (and inversion) OF VOLUMES
ROIOnly_MRI = ROIbox_MRI;
ROIOnly_PET = ROIbox_PET;
ROIOnly_MRI(maskBox_PET==0) = NaN;
ROIboxFill_MR = fillBox(ROIOnly_MRI);
ROIOnly_PET(maskBox_PET==0) = NaN;
minMRI = min(ROIOnly_MRI(:));
ROIboxFill_MR = ROIboxFill_MR - minMRI;
ROIOnly_MRI = ROIOnly_MRI - minMRI;
ROIboxFill_MR = ROIboxFill_MR./max(ROIOnly_MRI(:)).*255;
if strcmp(MRIinv,'Inv')
    ROIboxFill_MR = 255 - ROIboxFill_MR;
end
minPET = min(ROIOnly_PET(:));
ROIbox_PET = ROIbox_PET - minPET;
ROIOnly_PET = ROIOnly_PET - minPET;
ROIbox_PET = ROIbox_PET./max(ROIOnly_PET(:)).*255;


% WAVELET DECOMPOSITION AND FUSION OF COEFFICIENTS
wdecMRI = wavedec3(ROIboxFill_MR,1,wavelet);
wdecPET = wavedec3(ROIbox_PET,1,wavelet);
wdecFUSED = wdecPET;
nbcell = length(wdecMRI.dec);
for i = 1:nbcell
    wdecMRI.dec{i} = CollewetNorm(wdecMRI.dec{i});
    wdecMRI.dec{i}(isnan(wdecMRI.dec{i})) = wdecPET.dec{i}(isnan(wdecMRI.dec{i}));
    wdecFUSED.dec{i} = (1-MRIweight).*wdecPET.dec{i} + MRIweight.*wdecMRI.dec{i};
end


% WAVELET RECONSTRUCTION
fusedBox = waverec3(wdecFUSED);
fusedBox = fusedBox - min(fusedBox(:));
fusedBox = fusedBox./max(fusedBox(:)).*255;

end