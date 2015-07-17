function [sData] = readDICOMdir(dicomPath,waitB)
% -------------------------------------------------------------------------
% function [sData] = readDICOMdir(dicomPath,waitB)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function reads the DICOM content of a single directory. It then 
% organizes the data it in a cell of structures called 'sData', and 
% computes the region of interest (ROI) defined by a given RTstruct (if 
% present in the directory).
% -------------------------------------------------------------------------
% INPUTS:
% - dicomPath: Full path where the DICOM files to read are located.
% - waitB: Logical boolean. If true, a waiting bar will be displayed.
% -------------------------------------------------------------------------
% OUTPUTS:
% - sData: Cell of structures organizing the content of the volume data, 
%          DICOM headers, DICOM RTstruct* (used to compute the ROI) and 
%          DICOM REGstruct* (used to register a MRI volume to a PET volume)
%          * If present in the directory
%    --> sData{1}: Explanation of cell content
%    --> sData{2}: Imaging data and ROI defintion (if applicable)
%    --> sData{3}: DICOM headers of imaging data
%    --> sData{4}: DICOM RTstruct (if applicable)
%    --> sData{5}: DICOM REGstruct (if applicable)
% -------------------------------------------------------------------------
% AUTHOR(S): 
% - Martin Vallieres <mart.vallieres@gmail.com>
% - Sebastien Laberge <sebastien.laberge.3000@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: May 2015
%--------------------------------------------------------------------------
% STATEMENT:
% This file is part of <https://github.com/mvallieres/radiomics/>, 
% a package providing MATLAB programming tools for radiomics analysis.
% --> Copyright (C) 2015  Martin Vallieres, Sebastien Laberge
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


% INITIALIZATION
if waitB
    waitbarHandle = waitbar(0,'Loading DICOM files...','WindowStyle','modal');
end
elements = dir(dicomPath);
nElements = length(elements);
volume = cell(1,1,nElements);
dicomHeaders = [];
RTstruct = [];
REG = [];


% READING DIRECTORY CONTENT
sliceNumber = 0;
for elementNumber = 1:nElements
    elementName = elements(elementNumber).name;
    if ~strcmp(elementName,'.') && ~strcmp(elementName,'..') % Good enough for Linux, add conditions for MAC and Windows.
        elementFullFile = fullfile(dicomPath,elementName);
        if isdicom(elementFullFile)
            tmp = dicominfo(elementFullFile);
            if strcmp(tmp.Modality,'RTSTRUCT')
                RTstruct = tmp;
            elseif strcmp(tmp.Modality,'REG')
                REG = tmp;
            elseif strcmp(tmp.Modality,'MR') || strcmp(tmp.Modality,'PT') || strcmp(tmp.Modality,'CT')
                sliceNumber = sliceNumber + 1;
                volume{sliceNumber} = double(dicomread(elementFullFile));
                dicomHeaders = appendStruct(dicomHeaders,tmp);
            end
        end
    end
    if waitB
        waitbar(elementNumber/nElements,waitbarHandle); 
    end
end
nSlices = sliceNumber; % Total number of slices
volume = volume(1:nSlices); % Suppress empty cells in images


% DETERMINE THE SCAN ORIENTATION
dist = [abs(dicomHeaders(2).ImagePositionPatient(1) - dicomHeaders(1).ImagePositionPatient(1)), ...
        abs(dicomHeaders(2).ImagePositionPatient(2) - dicomHeaders(1).ImagePositionPatient(2)), ...
        abs(dicomHeaders(2).ImagePositionPatient(3) - dicomHeaders(1).ImagePositionPatient(3))];
[~,index] = max(dist);
if index == 1
    orientation = 'Sagittal';
elseif index == 2
    orientation = 'Coronal';
else
    orientation = 'Axial';
end


% SORT THE IMAGES AND DICOM HEADERS
slicePositions = zeros(1,nSlices);
for sliceNumber = 1:nSlices
    slicePositions(sliceNumber) = dicomHeaders(sliceNumber).ImagePositionPatient(index);
end
[~,indices] = sort(slicePositions);
volume = cell2mat(volume(indices));
dicomHeaders = dicomHeaders(indices);


% FILL sData
sData = cell(1,5);
type = dicomHeaders(1).Modality;
if strcmp(type,'PT') || strcmp(type,'CT')
    if strcmp(type,'PT')
        type = 'PET';
    end
    for i=1:size(volume,3)
        volume(:,:,i)=volume(:,:,i)*dicomHeaders(i).RescaleSlope + dicomHeaders(i).RescaleIntercept;
    end
end
type = [type,'scan'];
sData{1} = struct('Cell_1','Explanation of cell content', ...
                  'Cell_2','Imaging data and ROI defintion (if applicable)', ...
                  'Cell_3','DICOM headers of imaging data', ...
                  'Cell_4','DICOM RTstruct (if applicable)', ...
                  'Cell_5','DICOM REGstruct (if applicable)');

sData{2}.scan.volume = volume;
sData{2}.scan.orientation = orientation;
try sData{2}.scan.pixelW = dicomHeaders(1).PixelSpacing(1); catch sData{2}.scan.pixelW = []; end % Pixel Width
try sData{2}.scan.sliceT = dicomHeaders(1).SliceThickness; catch sData{2}.scan.sliceT = []; end % Slice Thickness
s1 = round(0.5*nSlices); s2 = round(0.5*nSlices) + 1; % Slices selected to calculate slice spacing
sData{2}.scan.sliceS = sqrt(sum((dicomHeaders(s1).ImagePositionPatient - dicomHeaders(s2).ImagePositionPatient).^2)); % Slice Spacing
sData{2}.type = type;
sData{3} = dicomHeaders;
sData{4} = RTstruct;
sData{5} = REG;


% COMPUTE TUMOR DELINEATION USING RTstruct
if ~isempty(sData{4})
    [sData] = computeROI(sData);
end


if waitB
    close(waitbarHandle)
end

end


% UTILITY FUNCTION
function [structureArray] = appendStruct(structureArray,newStructure)

if isempty(structureArray)
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

end