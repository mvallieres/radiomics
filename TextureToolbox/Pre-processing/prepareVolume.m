function [ROIonly,levels,ROIbox,maskBox] = prepareVolume(volume,mask,scanType,pixelW,sliceS,R,scale,textType,quantAlgo,Ng)
% -------------------------------------------------------------------------
% function [ROIonly,levels] = prepareVolume(volume,mask,scanType,pixelW,sliceS,R,scale,textType,quantAlgo,Ng)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function prepares the input volume for 3D texture analysis. The 
% following operations are performed:
%
% 1. Computation of the smallest box containing region of interest (ROI), 
%    if necessary (ROIbox).
% 2. Pre-processing of the ROIbox (PET: square-root, MR: Collewet 
%    normalizaton, CT: nothing).
% 3. Wavelet band-pass filtering (WBPF).
% 4. Isotropic resampling.
% 5. Quantization of intensity dynamic range.
%
% --> This function is compatible with 2D analysis (language not adapted in the text)
% -------------------------------------------------------------------------
% REFERENCE:
% [1] Vallieres, M. et al. (2015). A radiomics model from joint FDG-PET and 
%     MRI texture features for the prediction of lung metastases in soft-tissue 
%     sarcomas of the extremities. Physics in Medicine and Biology, 60(14), 
%     5471-5496. doi:10.1088/0031-9155/60/14/5471
% -------------------------------------------------------------------------
% INPUTS:
% - volume: 2D or 3D array containing the medical images to analyze
% - mask: 2D or 3D array of dimensions corresponding to 'volume'. The mask 
%         contains 1's in the region of interest (ROI), and 0's elsewhere.
% - scanType: String specifying the type of scan analyzed. Either 'PETscan', 
%             'MRscan' or 'Other'.
% - pixelW: Numerical value specifying the in-plane resolution (mm) of 'volume'.
% - sliceS: Numerical value specifying the slice spacing (mm) of 'volume'.
%           Put a random number for 2D analysis.
% - R: Numerical value specifying the ratio of weight to band-pass coefficients 
%      over the weigth of the rest of coefficients (HHH and LLL). Provide R=1 
%      to not perform wavelet band-pass filtering.    
% - scale: Numerical value specifying the scale at which 'volume' is isotropically 
%          resampled (mm). If a string 'pixelW' is entered as input, the
%          volume will be isotropically resampled at the initial in-plane
%          resolution of 'volume' specified by 'pixelW'.
% - textType: String specifying for which type of textures 'volume' is 
%             being prepared. Either 'Global' or 'Matrix'. If 'Global', the 
%             volume will be prepared for Global texture features computation. 
%             If 'Matrix',the volume will be prepared for matrix-based texture 
%             features computation (i.e. GLCM, GLRLM, GLSZM, NGTDM).
% - quantAlgo: String specifying the quantization algorithm to use on 'volume'. 
%              Either 'Equal' for equal-probability quantization, 'Lloyd'
%              for Lloyd-Max quantization, or 'Uniform' for uniform quantization.
%              Use only if textType is set to 'Matrix'.
% - Ng: Integer specifying the number of gray levels in the quantization process.
%       Use only if textType is set to 'Matrix'.
% -------------------------------------------------------------------------
% OUTPUTS:
% - ROIonly: Smallest box containing the ROI, with the imaging data of
%            the ready for texture analysis computations. Voxels outside
%            the ROI are set to NaNs.
% - levels: Vector containing the quantized gray-levels in the tumor region
%           (or reconstruction levels of quantization).
% - ROIbox: Smallest box containing the ROI. Optional output.
% - maskBox: Smallest mask containing the ROI. Optional output.
% -------------------------------------------------------------------------
% EXAMPLE:
% Let a PET scan be defined by 'volume', with 'mask' defining the ROI. The 
% PET scan has in-plane resolution of 4 mm, with slice spacing of 3.27 mm.
%             
% 1. To prepare 'volume' for matrix-based texture analysis at a scale of 
%    5 mm, without WBPF, using a Lloyd-Max quantization algorithm with 32 
%    gray-levels, run:
%
%    [ROIonly,levels] = prepareVolume(volume,mask,'PETscan',4,3.27,1,5,'Matrix','Lloyd',32)
% 
%    Next, use 'ROIonly' and 'levels' as inputs to 'getGLCM.m, 'getGLRLM.m',
%    'getGLSZM.m' or 'getNGTDM.m'.
%      
%
% 2. To prepare 'volume' for global texture analysis at a scale equal to the 
%    in-plane resolution, with R=2, run:
%
%    [ROIonly] = prepareVolume(volume,mask,'PETscan',4,3.27,2,'pixelW','Global')
% 
%    Next, use 'ROIonly' as input to 'getGlobalTextures.m'.
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


% VERIFICATION OF SOME INPUTS
if ~strcmp(textType,'Global') && ~strcmp(textType,'Matrix')
    error('Last argument must either be ''Global'' or ''Matrix''')
end
if nargin > 8
    if strcmp(quantAlgo,'Lloyd')
        quantization = @(x,y) lloydQuantization(x,y);
    elseif strcmp(quantAlgo,'Equal')
        quantization = @(x,y) equalQuantization(x,y);
    elseif strcmp(quantAlgo,'Uniform')
        quantization = @(x,y) uniformQuantization(x,y);
    else
        error('Error with quantization algorithm input. Must either be ''Equal'' or ''Lloyd'' or ''Uniform''')
    end
end


% COMPUTATION OF THE SMALLEST BOX CONTAINING THE ROI
[boxBound] = computeBoundingBox(mask);
maskBox = mask(boxBound(1,1):boxBound(1,2),boxBound(2,1):boxBound(2,2),boxBound(3,1):boxBound(3,2));
ROIbox = volume(boxBound(1,1):boxBound(1,2),boxBound(2,1):boxBound(2,2),boxBound(3,1):boxBound(3,2));


% PRE-PROCESSING OF ROI BOX
ROIbox = double(ROIbox);
if strcmp(scanType,'PETscan') || strcmp(scanType,'PTscan')
    ROIbox = sqrt(ROIbox);
elseif strcmp(scanType,'MRscan')
    ROIonly = ROIbox;
    ROIonly(~maskBox) = NaN;
    temp = CollewetNorm(ROIonly);
    maskBox(isnan(temp)) = 0;
end


% WAVELET BAND-PASS FILTERING
if R ~= 1
    if sum(isnan(ROIbox(:)))
        ROIbox = fillBox(ROIbox); % Necessary in cases we have a ROI box containing NaN's.
    end
    ROIbox = waveletBPfilt(ROIbox,R,'sym8');
end


% ISOTROPIC RESAMPLING
flagPW = 0;
if strcmp(scale,'pixelW')
    flagPW = 1;
end
if flagPW
    a = 1;
    b = 1;
    c = sliceS/pixelW;
else
    a = pixelW/scale;
    b = pixelW/scale;
    c = sliceS/scale;
end
if numel(size(ROIbox))==3
    if a + b + c ~= 3 % If false, no resampling is needed
        maskBox = imresize3D(maskBox,[],[round(double(size(maskBox,1))*a),round(double(size(maskBox,2))*b),round(double(size(maskBox,3))*c)],'nearest','fill');
        ROIbox = imresize3D(ROIbox,[],[round(double(size(ROIbox,1))*a),round(double(size(ROIbox,2))*b),round(double(size(ROIbox,3))*c)],'cubic','fill');
    end
elseif numel(size(ROIbox))==2
    if a + b ~= 2 % If false, no resampling is needed
        maskBox = imresize(maskBox,[round(double(size(maskBox,1))*a),round(double(size(maskBox,2))*b)],'nearest');
        ROIbox = imresize(ROIbox,[round(double(size(ROIbox,1))*a),round(double(size(ROIbox,2))*b)],'cubic','Antialiasing',true);
    end
end
ROIonly = ROIbox; ROIonly(~maskBox) = NaN; ROIonly(maskBox<0) = NaN;


% QUANTIZATION
if strcmp(textType,'Matrix')
    [ROIonly,levels] = quantization(ROIonly,Ng);
elseif strcmp(textType,'Global')
    levels = 0;
end

end
