function [textures] = getGLRLMtextures(GLRLM)
% -------------------------------------------------------------------------
% function [textures] = getGLRLMtextures(GLRLM)
% -------------------------------------------------------------------------
% DESCRIPTION:
% This function computes texture features from an input Gray-Level 
% Run-Length Matrix (GLRLM).
% -------------------------------------------------------------------------
% REFERENCES:
% [1] Galloway, M. M. (1975). Texture analysis using gray level run lengths. 
%     Computer Graphics and Image Processing, 4(2), 172–179.
% [2] Chu, A., Sehgal, C. M., & Greenleaf, J. F. (1990). Use of gray value 
%     distribution of run lengths for texture analysis. Pattern Recognition
%     Letters, 11(6), 415-419.
% [3] Dasarathy, B. V., & Holder, E. B. (1991). Image characterizations 
%     based on joint gray level-run length distributions. Pattern 
%     Recognition Letters, 12(8), 497-502.
% [4] Thibault, G., Fertil, B., Navarro, C., Pereira, S., Cau, P., Levy, 
%     N., Mari, J.-L. (2009). Texture Indexes and Gray Level Size Zone 
%     Matrix. Application to Cell Nuclei Classification. In Pattern 
%     Recognition and Information Processing (PRIP) (pp. 140–145).
% -------------------------------------------------------------------------
% INPUTS:
% - GLRLM: Gray-Level Run-Length Matrix.
%
% ** 'GLRLM' should be the output from 'getGLRLM.m' **
% -------------------------------------------------------------------------
% OUTPUTS:
% - textures: Struture specifying the values of different GLRLM texture
%             features as defined below.
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


% USEFUL MATRICES, VECTORS AND QUANTITIES
sz = size(GLRLM); % Size of GLRLM
nRuns = sum(GLRLM(:));
cVect = 1:sz(2); rVect = 1:sz(1);% Row and column vectors
[cMat,rMat] = meshgrid(cVect,rVect); % Column and row indicators for each entry of the GLRLM
pg = sum(GLRLM,2)'; % Gray-Level Run-Number Vector
pr = sum(GLRLM); % Run-Length Run-Number Vector


% COMPUTATION OF TEXTURE FEATURES
% 1. Short Run Emphasis (SRE), Ref.[1]
textures.SRE = (pr*(cVect.^(-2))')/nRuns;

% 2. Long Run Emphasis (LRE), Ref.[1]
textures.LRE = (pr*(cVect.^2)')/nRuns;

% 3. Gray-Level Nonuniformity (GLN), adapted from Ref.[1]
textures.GLN = sum(pg.^2)/nRuns;

% 4. Run-Length Nonuniformity (RLN), adapted from Ref.[1]
textures.RLN = sum(pr.^2)/nRuns;

% 5. Run Percentage (RP), adapted from Ref.[1]
textures.RP = nRuns/(pr*cVect');

% 6. Low Gray-Level Run Emphasis (LGRE), Ref.[2]
textures.LGRE = (pg*(rVect.^(-2))')/nRuns;

% 7. High Gray-Level Run Emphasis (HGRE), Ref.[2]
textures.HGRE = (pg*(rVect.^2)')/nRuns;

% 8. Short Run Low Gray-Level Emphasis (SRLGE), Ref.[3]
textures.SRLGE = sum(sum(GLRLM.*(rMat.^(-2)).*(cMat.^(-2))))/nRuns;

% 9. Short Run High Gray-Level Emphasis (SRHGE), Ref.[3]
textures.SRHGE = sum(sum(GLRLM.*(rMat.^2).*(cMat.^(-2))))/nRuns;

% 10. Long Run Low Gray-Level Emphasis (LRLGE), Ref.[3]
textures.LRLGE = sum(sum(GLRLM.*(rMat.^(-2)).*(cMat.^2)))/nRuns;

% 11. Long Run High Gray-Level Emphasis (LRHGE), Ref.[3]
textures.LRHGE = sum(sum(GLRLM.*(rMat.^2).*(cMat.^2)))/nRuns;


% New features according to Ref.[4]
GLRLM = GLRLM./nRuns;
pg = sum(GLRLM,2)'; pr = sum(GLRLM);
ug = (pg*rVect')/(sz(1)*sz(2));
ur = (pr*cVect')/(sz(1)*sz(2));

% 12. Gray-Level Variance (GLV), adapted from Ref.[4]
GLV = 0;
for g = 1:sz(1)
    for r = 1:sz(2)
        GLV = GLV + (GLRLM(g,r)*g-ug)^2;
    end
end
GLV = GLV/(sz(1)*sz(2));
textures.GLV = GLV;

% 13. Run-Length Variance (RLV), adapted from Ref.[4]
RLV = 0;
for g = 1:sz(1)
    for r = 1:sz(2)
        RLV = RLV + (GLRLM(g,r)*r-ug)^2;
    end
end
RLV = RLV/(sz(1)*sz(2));
textures.RLV = RLV;

end