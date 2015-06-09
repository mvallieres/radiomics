function [testSets] = findBootTestSet(bootSam)
% -------------------------------------------------------------------------
% function [testSet] = findBootTestSet(trainSet)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function finds the testing sets for a set of bootstrap samples 
% (i.e. of training sets). Testing sets for different bootstrap samples
% may not be of the same length, so the output is a cell of testing sets.
% -------------------------------------------------------------------------
% INPUTS:                             
% - bootSam: Matrix of size [nInst X nBoot], where nInst is the number of
%            instances in the outcome vector.
% -------------------------------------------------------------------------
% OUTPUTS:
% - testSets: Cell of size [1 X nBoot], where each entry is a testing set
%             column vector.
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

nInst = size(bootSam,1);
nBoot = size(bootSam,2);
testSets = cell(1,nBoot);
for n = 1:nBoot
    vectTest = ones(nInst,1);
    vectTest(bootSam(:,n)) = 0;
    testSets{n} = find(vectTest);
end

end