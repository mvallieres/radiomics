function [XtrainOut,XtestOut] = normalizeZeroOne(XtrainIn,XtestIn)
% -------------------------------------------------------------------------
% function [XtrainOut,XtestOut] = normalizeZeroOne(XtrainIn,XtestIn)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function normalizes each feature of an input training matrix between
% 0 and 1, and optionally, applies the same normalization to an input 
% testing matrix.
% -------------------------------------------------------------------------
% INPUTS:                             
% - XtrainIn: Matrix of size [nInst X nFeat], specifying the numerical data 
%             of the features of the input features, where 'nInst' refers 
%             to the number of instances in XtrainIn, and 'nFeat' to the 
%             number of features in XtrainIn. Each column is a different 
%             feature.
% - XtestIn: (optional input). Matrix of size [nInst X nFeat], specifying 
%             the numerical data of the features of the input testing data. 
%             Each column is a different feature.
% -------------------------------------------------------------------------
% OUTPUTS:
% - XtrainOut: Normalized XtrainIn.
% - XtestOut: (optional output). Normalized XtestIn.
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


szTrain = size(XtrainIn);
XtrainOut = XtrainIn;
if nargout == 2 && nargin == 2
    XtestOut = XtestIn;
end

for i = 1:szTrain(2)
    minVal = min(XtrainIn(:,i));
    XtrainOut(:,i) = XtrainOut(:,i)-minVal;
    maxVal = max(XtrainOut(:,i));
    XtrainOut(:,i) = XtrainOut(:,i)./maxVal;
    if nargout == 2
        XtestOut(:,i)=XtestOut(:,i)-minVal;
        XtestOut(:,i)=XtestOut(:,i)./maxVal;
    end
end

end