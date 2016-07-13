function [textures] = getGLCMtextures(GLCM)
% -------------------------------------------------------------------------
% function [textures] = getGLCMtextures(GLCM))
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function computes texture features from an input Gray-Level 
% Co-occurence Matrix (GLCM).
% -------------------------------------------------------------------------
% REFERENCES:
% [1] Haralick, R. M., Shanmugam, K., & Dinstein, I. (1973). Textural 
%     features for image classification. IEEE Transactions on Systems, Man 
%     and Cybernetics, smc 3(6), 610–621.
% [2] Assefa, D., Keller, H., Ménard, C., Laperriere, N., Ferrari, R. J., & 
%     Yeung, I. (2010). Robust texture features for response monitoring of 
%     glioblastoma multiforme on T1 -weighted and T2 -FLAIR MR images: A 
%     preliminary investigation in terms of identification and segmentation. 
%     Medical Physics, 37(4), 1722–1736.
% [3] Thibault, G. (2009). Indices de formes et de textures: de la 2D vers 
%     la 3D. Application au classement de noyaux de cellules. PhD Thesis, 
%     Université AIX-Marseille: p.172.
% [4] Aerts, H.J.W.L. et al. Decoding tumour phenotype by noninvasive
%     imaging using a quantitative radiomics approach. Nat. Commun. 5:4006
%     doi: 10.1038/ncomms5006 (2014).
% -------------------------------------------------------------------------
% INPUTS:
% - GLCM: Gray-Level Co-occurence Matrix.
%
% ** 'GLCM' should be the output from 'getGLCM.m' **
% -------------------------------------------------------------------------
% OUTPUTS:
% - textures: Struture specifying the values of different GLCM texture
%             features as defined below.
% -------------------------------------------------------------------------
% AUTHOR(S): Martin Vallieres, <mart.vallieres@gmail.com>
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
textures = struct;
matrixtemp = GLCM;
GLCM = GLCM/(sum(GLCM(:))); % Normalization of GLCM
nL = max(size(GLCM));
indVect = 1:nL;
[colGrid,rowGrid] = meshgrid(indVect,indVect);


% COMPUTATION OF TEXTURE FEATURES
% 1. Energy, Ref.[1]
textures.Energy = sum(sum(GLCM.^2));

% 2. Contrast, Ref.[1]
contrast = 0.0;
for n = 0:nL-1
   temp = 0;
   for i = 1:nL
      for j = 1:nL
         if (abs(i-j) == n)
            temp = temp+GLCM(i,j);
         end
      end
   end
   contrast = contrast + n^2*temp;
end
textures.Contrast = contrast;

% 3. Entropy, Ref.[1]
textures.Entropy = -sum(sum(GLCM.*log2(GLCM + realmin)));

% 4. Homogeneity, adapted from Ref.[1]
temp = 0;
for i = 1:nL
   for j = 1:nL
      temp = temp + GLCM(i,j)/(1+abs(i-j));
   end
end
textures.Homogeneity = temp;

% 5. Correlation, adapted from Ref. [1] (this definition from MATLAB is preferred from the original one in [1])
textures.Correlation = graycoprops(round(matrixtemp),'Correlation');
textures.Correlation = struct2cell(textures.Correlation);
textures.Correlation = textures.Correlation{1};

% 6. Variance, Ref.[2]; and 7. SumAverage, Ref.[2]. (adapted from Variance and SumAverage metrics defined by Haralick in Ref. [1])
% However, in order to compare GLCMs of different sizes, the metrics
% are divided by the total number of elements in the GLCM (nL*nL). Also,
% there is probably an error in Assefa's paper [2]: in the variance equation,
% 'u' should appropriately be replaced by 'ux' and 'uy' as calculated in A1
% and A2 of the same paper (read 'ui' and 'uj' in our code).
ui = indVect*sum(GLCM,2);
uj = indVect*sum(GLCM)';
tempS = rowGrid.*GLCM + colGrid.*GLCM;
tempV = (rowGrid-ui).^2.*GLCM + (colGrid-uj).^2.*GLCM;
textures.SumAverage = 0.5*sum(tempS(:))/(nL^2);
textures.Variance = 0.5*sum(tempV(:))/(nL^2);

% 8. Dissimilarity, Ref.[3] 
diffMat = abs(rowGrid-colGrid);
temp = diffMat.*GLCM;
textures.Dissimilarity = sum(temp(:));

% 9. AutoCorrelation, Ref.[4]
temp = rowGrid .* colGrid .* GLCM;
textures.AutoCorrelation = sum(temp(:));


end