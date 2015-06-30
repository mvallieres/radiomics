function [textCellOut] = removeParam(textCellIn,paramUsed,baseline)
% -------------------------------------------------------------------------
% function [textCellOut] = removeParam(textCellIn,paramUsed,baseline)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function removes textures for which the extraction parameters are 
% not allowed to vary as defined by 'paramUsed', and only keep baseline 
% extraction parameters as defined by 'baseline'. See ref. [1] for more 
% details.
% -------------------------------------------------------------------------
% REFERENCE:
% [1] Vallieres, M. et al. (2015). A radiomics model from joint FDG-PET and 
%     MRI texture features for the prediction of lung metastases in soft-tissue 
%     sarcomas of the extremities. Physics in Medicine and Biology, 60(14), 
%     5471-5496. doi:10.1088/0031-9155/60/14/5471 
% -------------------------------------------------------------------------
% INPUTS:                             
% - textCellIn: Format of the files saved as output to organizeSeparateTextures_STS.m 
%               or organizeFusedTextures_STS.m, for example.
% - paramUsed: Vector of 1's and 0's to specify the degree of freedom on 
%              texture extraction parameters for the current experiment. 
%              For example, for an experiment where extraction parameters 
%              1, 2 and 4 in paramAll are allowed to vary, use
%              paramUsed = [1,1,0,1]. The maximal length of paraUsed
%              currently supported is 6.
% - baseline: Vector of numerical values specifying the baseline texture 
%             extraction parameters for each entry in paramUsed. See EXAMPLE
%             below for more details.
%
% See <https://github.com/mvallieres/radiomics/tree/master/STUDIES/STSstudy/Functions>
% to find organizeSeparateTextures_STS.m organizeFusedTextures_STS.m. 
% See masterScript_STS.m for a complete example of how to use the current 
% function.
% -------------------------------------------------------------------------
% OUTPUTS:
% - textCellOut: Modified textCellIn.
% -------------------------------------------------------------------------
% EXAMPLE:
% MRIinv_cell = {'NoInv','Inv'};
% MRIweight_mat = [1/4,1/3,1/2,2/3,3/4];
% R_mat = [1/2,2/3,1,3/2,2];
% scale_cell = {'pixelW',1,2,3,4,5};
% algo_cell = {'Equal','Lloyd'};
% Ng_mat = [8,16,32,64];
% paramAll = {MRIinv_cell,MRIweight_mat,R_mat,scale_cell,algo_cell,Ng_mat};
% paramUsed = [1 1 0 0 0 1];
% baseline = [1 3 3 1 2 3];
% -------------------------------------------------------------------------
% AUTHOR(S): Martin Vallieres <mart.vallieres@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: May 2015
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

textCellOut = textCellIn;
vectOut = find(~paramUsed);
nOut = length(vectOut);
for i = 1:nOut
    keep = baseline(vectOut(i));
    if vectOut(i) == 1
        for j = 1:keep-1
            textCellOut(1,:,:,:,:,:) = [];
        end
        for j = (keep+1):size(textCellIn,vectOut(i))
            textCellOut(end,:,:,:,:,:) = [];
        end
    elseif vectOut(i) == 2
        for j = 1:keep-1
            textCellOut(:,1,:,:,:,:) = [];
        end
        for j = (keep+1):size(textCellIn,vectOut(i))
            textCellOut(:,end,:,:,:,:) = [];
        end
    elseif vectOut(i) == 3
        for j = 1:keep-1
            textCellOut(:,:,1,:,:,:) = [];
        end
        for j = (keep+1):size(textCellIn,vectOut(i))
            textCellOut(:,:,end,:,:,:) = [];
        end
    elseif vectOut(i) == 4
        for j = 1:keep-1
            textCellOut(:,:,:,1,:,:) = [];
        end
        for j = (keep+1):size(textCellIn,vectOut(i))
            textCellOut(:,:,:,end,:,:) = [];
        end
    elseif vectOut(i) == 5
        for j = 1:keep-1
            textCellOut(:,:,:,:,1,:) = [];
        end
        for j = (keep+1):size(textCellIn,vectOut(i))
            textCellOut(:,:,:,:,end,:) = [];
        end
    elseif vectOut(i) == 6
        for j = 1:keep-1
            textCellOut(:,:,:,:,:,1) = [];
        end
        for j = (keep+1):size(textCellIn,vectOut(i))
            textCellOut(:,:,:,:,:,end) = [];
        end
    end
end

end