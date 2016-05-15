function [Xout,Yout] = shufflePartition(Xin,Yin)
% -------------------------------------------------------------------------
% function [Xout,Yout] = shufflePartition(Xin,Yin)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function shuffles the rows of a matrix 'Xin' and a column vector
% 'Yin' (same permutations to both variables).
% -------------------------------------------------------------------------
% INPUTS:                             
% - Xin: Input 2D matrix
% - Yin: Input column vector (optional)
% -------------------------------------------------------------------------
% OUTPUTS:
% - Xin: Output shuffled 2D matrix (shuffled row-wise)
% - Yin: Output shuffled column vector (optional)
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

% RANDOM NUMBER GENERATOR SEED
if ~RandStream.getGlobalStream.Seed
    rng('shuffle')
end

n = size(Xin,1);
nF = size(Xin,2);
if isa(Xin,'table')
    varName = Xin.Properties.VariableNames;
    Xin = table2cell(Xin);
    Xout = cell(n,nF);
    table = true;
else
    Xout = zeros(n,nF);
    table = false;
end

if nargout == 2
    Yout = zeros(n,1);
end

for i = 1:n
    index = ceil((n+1-i)*rand(1,1));
    Xout(i,1:nF) = Xin(index,1:nF);
    Xin(index,:) = [];
    if nargout == 2
        Yout(i,1) = Yin(index,1);
    end
    if nargin == 2
        Yin(index) = [];
    end
end

if table
    Xout = cell2table(Xout,'VariableNames',varName);
end

end