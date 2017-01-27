function [prob] = oobPredict_table(RF,X)
% -------------------------------------------------------------------------
% function [prob] = oobPredict_table(RF,X)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function computes out-of-bag predictions.
% -------------------------------------------------------------------------
% INPUTS:
% 1. RF: Random forest object out of TreeBagger.m
% 2. X: Array of size [nPatient X nFeatures] in "Table" format.
%       --> Ex: See dataStruct.mat created in organizeRF_LGG.m
% -------------------------------------------------------------------------
% OUTPUTS: 
% 1. prob: Vector of probability of prediction of an outcome for 
%          all patients.
% -------------------------------------------------------------------------
% AUTHOR(S): Martin Vallieres <mart.vallieres@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: January 2017
%--------------------------------------------------------------------------
% STATEMENT:
% This file is part of <https://github.com/mvallieres/radiomics/>, 
% a package providing MATLAB programming tools for radiomics analysis.
% --> Copyright (C) 2015-2017  Martin Vallieres
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
prob = zeros(nInst,1);
oobIndAll = RF.OOBIndices;
for i = 1:nInst
    oobInd = find(oobIndAll(i,:)); nTrees = numel(oobInd);
    pred = 0;
    for t = 1:nTrees
        pred = pred + str2double(RF.Trees{oobInd(t)}.predict(X(i,:)));
    end
    prob(i) = pred/nTrees;
end

end