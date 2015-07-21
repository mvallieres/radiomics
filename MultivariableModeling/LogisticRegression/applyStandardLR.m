function [coeff] = applyStandardLR(X,Y)
% -------------------------------------------------------------------------
% function [coeff] = applyStandardLR(X,Y)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function computes standard logistic regression using the utilities 
% from DREES <http://www.cerr.info/drees>.
% -------------------------------------------------------------------------
% INPUTS:
% - X: Matrix of size [nInst X nFeat], specifying the numerical data of the 
%      features of the input features, where 'nInst' refers to the number 
%      of instances in X, and 'nFeat' to the number of features in X. 
%      Each column is a different feature.
% - Y: Column vector of size [nInst X 1] specifying the outcome status 
%      (1 or 0) for all instances.
% -------------------------------------------------------------------------
% OUTPUTS:
% - coeff: Column vector of size [nCoeff+1 X 1] specifying the final 
%          logistic regression coefficients. Last entry specify the offset
%          of the model.
% -------------------------------------------------------------------------
% AUTHOR(S): 
% - DREES development team <http://www.cerr.info/drees> (logistic regression)
% - Martin Vallieres <mart.vallieres@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: July 2015
%--------------------------------------------------------------------------
% STATEMENT:
% This file is part of <https://github.com/mvallieres/radiomics/>, 
% a package providing MATLAB programming tools for radiomics analysis.
% --> Copyright 2010, Joseph O. Deasy, on behalf of the DREES development team.
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

 [~,coeff,~] = drxlr_apply_logistic_regression(X,Y);

end