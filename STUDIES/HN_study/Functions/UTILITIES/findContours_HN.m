function [contourVect] = findContours_HN(sData,cellName,scanType)
% -------------------------------------------------------------------------
% function [contourVect] = findContours_HN(sData,cellName)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% Finds the contour numbers in the sData structure according to the
% 'cellName' variable.
% -------------------------------------------------------------------------
% INPUTS:
% 1. sData: Cell of structures organizing the data.
% 2. cellName: Cell of strings of size {1,nConcatenate} specifying all the 
%              contour numbers to find in sData.
%              --> Ex: {'GTV1,GTV2','GTVn'}. For each entry, the function 
%              will find all associated contour numbers in sData. In the 
%              example above, the function will find 3 contour numbers 
%              (comma + cell1/cell2 separation). Different contours in a
%              single entry must be separated by one comma, and nothing
%              else. 'nConcatenate' can range from 1 to infinity.
% 3. scanType: (optional). String specifying the type of scan the contour 
%              belongs to. Either 'CT' or 'PET'. Default is 'CT'.
%              --> Ex: 'PET'
% -------------------------------------------------------------------------
% OUTPUTS:
% 1. contourVect: Numerical array of size [1,nContour] specifying all the 
%                 contour numbers found in sData.
%                 --> Ex: [1,4,3]
% -------------------------------------------------------------------------
% AUTHOR(S): 
% - Martin Vallieres <mart.vallieres@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: March 2016
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


contourNames = {};
nEntry = numel(cellName);
for i = 1:nEntry
    string = cellName{i};
    indComma = strfind(string,',');
    if isempty(indComma)
        contourNames = [contourNames,string];
    else
        indComma = [0,indComma,numel(string)+1];
        for j = 1:numel(indComma)-1
            contourNames = [contourNames;string((indComma(j)+1):((indComma(j+1)-1)))];
        end
    end
end

nContour = numel(contourNames);
contourVect = zeros(1,nContour);
nContourData = numel(sData{2}.scan.contour);
for i = 1:nContour
    for j = 1:nContourData
        name = sData{2}.scan.contour(j).name;
        if strcmp(scanType,'PET')
            if strcmp(name,['RTstruct_CTsim->PET(PET-CT)--',contourNames{i}])
                contourVect(i) = j;
                break
            end
        elseif strcmp(scanType,'CT')
            if strcmp(name,['RTstruct_CTsim->CT(PET-CT)--',contourNames{i}])
                contourVect(i) = j;
                break
            elseif ~isempty(strfind(name,contourNames{i})) && isempty(strfind(name,'RTstruct_CTsim->PET(PET-CT)--')) % In case the wanted contour (RTstruct_CTsim->CT(PET-CT)) is not found, better to directly use the contour from RT planning, if present.
                ind = strfind(name,'--');
                name = name((ind+2):end);
                if strcmp(name,contourNames{i})
                    contourVect(i) = j;
                end
            end
        end
    end
end

end