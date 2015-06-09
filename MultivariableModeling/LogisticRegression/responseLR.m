function [resp] = responseLR(X,coeff)
% -------------------------------------------------------------------------
% function [resp] = responseLR(X,coeff)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function computes the multivariable response of an input set of
% features from a given set of logistic regression coefficients.
% -------------------------------------------------------------------------
% INPUTS:                             
% - X: Matrix of size [nInst X nFeat], specifying the numerical data of the
%      features of the input training data, where 'nInst' refers to the 
%      number of instances in X, and 'nFeat' to the number of features in 
%      X. Each column is a different feature.
% - coeff: Column vector of size [nFeat+1 X 1] representing the set of 
%          logistic regression coefficients computed from
%          drxlr_apply_logistic_regression.m. One coefficient is present
%          for each feature in 'X', in addition to one offset coefficient.
% -------------------------------------------------------------------------
% OUTPUTS:
% - resp: Multivariable response vector of size [nInst X 1].
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

nInst = size(X,1);
nFeat = size(X,2);
resp = zeros(nInst,1);
for j = 1:nFeat
    resp(:) = resp(:) + X(:,j).*coeff(j);
end
resp(:) = resp(:) + coeff(end);

end