function [GLRLM] = getGLRLM(ROIonly,levels)
% -------------------------------------------------------------------------
% function [GLRLM] = getGLRLM(ROIonly,levels)
% -------------------------------------------------------------------------
% DESCRIPTION:
% This function computes the Gray-Level Run-Length Matrix (GLRLM) of the 
% region of interest (ROI) of an input volume. The input volume is assumed 
% to be isotropically resampled. Only one GLRLM is computed per scan, 
% simultaneously adding up all possible run-lengths in the 13 directions of 
% the 3D space. To account for discretization length differences, runs 
% constructed from voxels separated by a distance of sqrt(3) increment the 
% GLRLM by a value of sqrt(3), runs constructed from voxels separated by a 
% distance of sqrt(2) increment the GLRLM by a value of sqrt(2), and runs 
% constructed from voxels separated by a distance of 1 increment the GLRLM 
% by a value of 1. This function uses other functions from Wei's GLRLM 
% toolbox [2].
%
% --> This function is compatible with 2D analysis (language not adapted in the text)
% -------------------------------------------------------------------------
% REFERENCES:
% [1] Galloway, M. M. (1975). Texture analysis using gray level run lengths.
%     Computer Graphics and Image Processing, 4(2), 172–179.
% [2] Wei's GLRLM toolbox: Xunkai Wei, Gray Level Run Length Matrix Toolbox
%     v1.0, Software,Beijing Aeronautical Technology Research Center, 2007.
%     <http://www.mathworks.com/matlabcentral/fileexchange/17482-gray-level-run-length-matrix-toolbox>
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
% - GLRLM: Gray-Level Run-Length Matrix of 'ROIOnly'.
% -------------------------------------------------------------------------
% AUTHOR(S): 
% - Martin Vallieres <mart.vallieres@gmail.com>
% - Xunkai Wei <xunkai.wei@gmail.com>
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
%
%    _______________________________________________________________
%
% --> Copyright (c) 2007-2012, Xunkai Wei
%     All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%      
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.
% -------------------------------------------------------------------------


% PRELIMINARY
nLevel = length(levels);
if nLevel > 100
    adjust = 10000;
else
    adjust = 1000;
end
levelTemp = max(levels)+1;
ROIonly(isnan(ROIonly)) = levelTemp; % Last row needs to be taken out of the GLRLM
levels = [levels,levelTemp];


% QUANTIZATION EFFECTS CORRECTION
% In case (for example) we initially wanted to have 64 levels, but due to
% quantization, only 60 resulted.
uniqueVol = round(levels*adjust)/adjust;
ROIonly=round(ROIonly*adjust)/adjust;
NL = length(levels) - 1;


%INITIALIZATION
sizeV = size(ROIonly);
numInit = ceil(max(sizeV)*sqrt(3)); % Max run length
GLRLM = zeros(NL+1,numInit);


% START COMPUTATION
% Directions [1,0,0], [0 1 0], [1 1 0] and [-1 1 0] : 2D directions
% (x:right-left, y:top-bottom, z:3rd dimension)  
if numel(size(ROIonly)) == 3
    nComp = sizeV(3); % We can add-up the GLRLMs taken separately in every image in the x-y plane
else
    nComp = 1;
end
for i = 1:nComp
    image = ROIonly(:,:,i);
    uniqueIm = unique(image);
    NLtemp = length(uniqueIm);
    indexRow = zeros(NLtemp,1);
    temp = image;
    for j = 1:NLtemp
        indexRow(j) = find(uniqueIm(j)==uniqueVol);
        image(temp==uniqueIm(j)) = j;
    end
    
    % [1,0,0]
    GLRLMtemp = rle_0(image,NLtemp);
    nRun = size(GLRLMtemp,2);
    GLRLM(indexRow(1:NLtemp),1:nRun) = GLRLM(indexRow(1:NLtemp),1:nRun) + GLRLMtemp(1:NLtemp,1:nRun); % Cumulative addition into the GLRLM
    
    % [0 1 0]
    GLRLMtemp = rle_0(image',NLtemp);
    nRun = size(GLRLMtemp,2);
    GLRLM(indexRow(1:NLtemp),1:nRun) = GLRLM(indexRow(1:NLtemp),1:nRun) + GLRLMtemp(1:NLtemp,1:nRun); % Cumulative addition into the GLRLM
    
    % [1 1 0]
    seq = zigzag(image);
    GLRLMtemp = rle_45(seq,NLtemp);
    nRun = size(GLRLMtemp,2);
    GLRLM(indexRow(1:NLtemp),1:nRun) = GLRLM(indexRow(1:NLtemp),1:nRun) + GLRLMtemp(1:NLtemp,1:nRun)*sqrt(2); % Cumulative addition into the GLRLM
    
    % [-1 1 0]
    seq = zigzag(fliplr(image));
    GLRLMtemp = rle_45(seq,NLtemp);
    nRun = size(GLRLMtemp,2);
    GLRLM(indexRow(1:NLtemp),1:nRun) = GLRLM(indexRow(1:NLtemp),1:nRun) + GLRLMtemp(1:NLtemp,1:nRun)*sqrt(2); % Cumulative addition into the GLRLM
end

if numel(size(ROIonly)) == 3 % 3D DIRECTIONS
    % Directions [0,0,1], [1 0 1] and [-1 0 1]
    % (x:right-left, y:top-bottom, z:3rd dimension)
    nComp = sizeV(1); % We can add-up the GLRLMs taken separately in every image in the x-z plane
    image = zeros(sizeV(3),sizeV(2));
    for i = 1:nComp
        for j = 1:sizeV(3)
            image(j,1:end) = ROIonly(i,1:end,j);
        end
        uniqueIm = unique(image);
        NLtemp = length(uniqueIm);
        indexRow = zeros(NLtemp,1);
        temp = image;
        for j=1:NLtemp
            indexRow(j) = find(uniqueIm(j)==uniqueVol);
            image(temp==uniqueIm(j)) = j;
        end
        
        % [0,0,1]
        GLRLMtemp = rle_0(image',NLtemp);
        nRun = size(GLRLMtemp,2);
        GLRLM(indexRow(1:NLtemp),1:nRun) = GLRLM(indexRow(1:NLtemp),1:nRun) + GLRLMtemp(1:NLtemp,1:nRun); % Cumulative addition into the GLRLM
        
        % [1 0 1]
        seq = zigzag(image);
        GLRLMtemp = rle_45(seq,NLtemp);
        nRun = size(GLRLMtemp,2);
        GLRLM(indexRow(1:NLtemp),1:nRun) = GLRLM(indexRow(1:NLtemp),1:nRun) + GLRLMtemp(1:NLtemp,1:nRun)*sqrt(2); % Cumulative addition into the GLRLM
        
        % [-1 0 1]
        seq = zigzag(fliplr(image));
        GLRLMtemp = rle_45(seq,NLtemp);
        nRun = size(GLRLMtemp,2);
        GLRLM(indexRow(1:NLtemp),1:nRun) = GLRLM(indexRow(1:NLtemp),1:nRun) + GLRLMtemp(1:NLtemp,1:nRun)*sqrt(2); % Cumulative addition into the GLRLM
    end

    % Directions [0,1,1] and [0 -1 1]
    % (x:right-left, y:top-bottom, z:3rd dimension)
    nComp = sizeV(2); % We can add-up the GLRLMs taken separately in every image in the y-z plane
    image = zeros(sizeV(1),sizeV(3));
    for i = 1:nComp
        for j = 1:sizeV(3)
            image(1:end,j) = ROIonly(1:end,i,j);
        end
        uniqueIm = unique(image);
        NLtemp = length(uniqueIm);
        indexRow = zeros(NLtemp,1);
        temp = image;
        for j = 1:NLtemp
            indexRow(j) = find(uniqueIm(j)==uniqueVol);
            image(temp==uniqueIm(j)) = j;
        end
        
        % [0,1,1]
        seq = zigzag(image);
        GLRLMtemp = rle_45(seq,NLtemp);
        nRun = size(GLRLMtemp,2);
        GLRLM(indexRow(1:NLtemp),1:nRun) = GLRLM(indexRow(1:NLtemp),1:nRun) + GLRLMtemp(1:NLtemp,1:nRun)*sqrt(2); % Cumulative addition into the GLRLM
        
        % [0 -1 1]
        seq = zigzag(fliplr(image));
        GLRLMtemp = rle_45(seq,NLtemp);
        nRun = size(GLRLMtemp,2);
        GLRLM(indexRow(1:NLtemp),1:nRun) = GLRLM(indexRow(1:NLtemp),1:nRun) + GLRLMtemp(1:NLtemp,1:nRun)*sqrt(2); % Cumulative addition into the GLRLM
    end

    % Four corners: [1,1,1], [-1,1,1], [-1,1,-1], [1,1,-1]
    % (x:right-left, y:top-bottom, z:3rd dimension)
    image = zeros(sizeV(3),sizeV(2));
    temp = rand(sizeV(3),sizeV(2));
    diagTemp = spdiags(temp);
    szDiag = size(diagTemp);
    diagMat1 = zeros(szDiag(1),szDiag(2),sizeV(1));
    diagMat2 = zeros(size(diagTemp,1),size(diagTemp,2),sizeV(1));
    for i = 1:sizeV(1)
        for j = 1:sizeV(3)
            image(j,1:end) = ROIonly(i,1:end,j);
        end
        try
            diagMat1(:,:,i)=spdiags(image);
        catch
            % Add a column at the beginning to prevent errors
            temp=spdiags(image);
            numberDiff=abs(size(temp,2)-size(diagMat1,2));
            if mod(numberDiff,2) % Odd difference number
                temp=padarray(temp,[0,(numberDiff+1)/2,0],0);
                diagMat1(:,:,i)=temp(:,1:end-1);
            else
                diagMat1(:,:,i)=padarray(temp,[0,numberDiff/2,0],0);
            end
        end
        try
            diagMat2(:,:,i)=spdiags(fliplr(image));
        catch
            % Add a column at the beginning to prevent errors
            temp = spdiags(fliplr(image));
            numberDiff = abs(size(temp,2)-size(diagMat2,2));
            if mod(numberDiff,2) % Odd difference number
                temp = padarray(temp,[0,(numberDiff+1)/2,0],0);
                diagMat2(:,:,i) = temp(:,1:end-1);
            else
                diagMat2(:,:,i) = padarray(temp,[0,numberDiff/2,0],0);
            end
        end
    end
    for j = 1:szDiag(2)
        index = (diagMat1(:,j,1)~=0);
        nTemp = sum(index);
        image1 = zeros(sizeV(1),nTemp);
        image2 = zeros(sizeV(1),nTemp);
        for k = 1:sizeV(1)
            image1(k,1:nTemp) = diagMat1(index(1:end),j,k)';
            image2(k,1:nTemp) = diagMat2(index(1:end),j,k)';
        end
        
        % 2 first corners
        uniqueIm = unique(image1);
        NLtemp = length(uniqueIm);
        indexRow = zeros(NLtemp,1);
        temp = image1;
        for i = 1:NLtemp
            indexRow(i) = find(uniqueIm(i)==uniqueVol);
            image1(temp==uniqueIm(i)) = i;
        end
        seq = zigzag(image1);
        GLRLMtemp = rle_45(seq,NLtemp);
        nRun = size(GLRLMtemp,2);
        GLRLM(indexRow(1:NLtemp),1:nRun) = GLRLM(indexRow(1:NLtemp),1:nRun) + GLRLMtemp(1:NLtemp,1:nRun)*sqrt(3); % Cumulative addition into the GLRLM
        seq = zigzag(fliplr(image1));
        GLRLMtemp = rle_45(seq,NLtemp);
        nRun = size(GLRLMtemp,2);
        GLRLM(indexRow(1:NLtemp),1:nRun) = GLRLM(indexRow(1:NLtemp),1:nRun) + GLRLMtemp(1:NLtemp,1:nRun)*sqrt(3); % Cumulative addition into the GLRLM
        
        % 2 last corners
        uniqueIm = unique(image2);
        NLtemp = length(uniqueIm);
        indexRow = zeros(NLtemp,1);
        temp = image2;
        for i = 1:NLtemp
            indexRow(i) = find(uniqueIm(i)==uniqueVol);
            image2(temp==uniqueIm(i)) = i;
        end
        seq = zigzag(image2);
        GLRLMtemp = rle_45(seq,NLtemp);
        nRun = size(GLRLMtemp,2);
        GLRLM(indexRow(1:NLtemp),1:nRun) = GLRLM(indexRow(1:NLtemp),1:nRun) + GLRLMtemp(1:NLtemp,1:nRun)*sqrt(3); % Cumulative addition into the GLRLM
        seq = zigzag(fliplr(image2));
        GLRLMtemp = rle_45(seq,NLtemp);
        nRun = size(GLRLMtemp,2);
        GLRLM(indexRow(1:NLtemp),1:nRun) = GLRLM(indexRow(1:NLtemp),1:nRun) + GLRLMtemp(1:NLtemp,1:nRun)*sqrt(3); % Cumulative addition into the GLRLM
    end
end


% REMOVE UNECESSARY COLUMNS
GLRLM(end,:) = [];
stop = find(sum(GLRLM),1,'last');
GLRLM(:,(stop+1):end) = [];

end