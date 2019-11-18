function [structureArray] = appendStruct(structureArray,newStructure)
% -------------------------------------------------------------------------
% function [structureArray] = appendStruct(structureArray,newStructure)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function apprends a structure to an already existing one.
% -------------------------------------------------------------------------
% INPUTS:
% - structureArray: Initial structure.
% - newStructure: New structure to append to structureArray.
% -------------------------------------------------------------------------
% OUTPUTS:
% - structureArray: Appended structure array.
% -------------------------------------------------------------------------
% AUTHOR(S): 
% - Sebastien Laberge <sebastien.laberge.3000@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: January 2016
%--------------------------------------------------------------------------
% STATEMENT:
% This file is part of <https://github.com/mvallieres/radiomics/>, 
% a package providing MATLAB programming tools for radiomics analysis.
% --> Copyright (C) 2015  Sebastien Laberge
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

if isempty(structureArray) || isempty(fieldnames(structureArray))
    structureArray = newStructure;
    return
end

structLength = length(structureArray);
fields = fieldnames(structureArray(1));
nFields = length(fields);

for i = 1:nFields
    try
        structureArray(structLength + 1).(fields{i}) = newStructure.(fields{i});
    catch
        structureArray(structLength + 1).(fields{i}) = 'FIELD NOT PRESENT';
    end
end