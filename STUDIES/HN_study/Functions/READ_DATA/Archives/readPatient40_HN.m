function [sDataPET,sDataCT] = readPatient40_HN(pathPET,pathCT)
% -------------------------------------------------------------------------
% function [sDataPET,sDataCT] = readPatient40_HN(pathPET,pathCT)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% Reason of existence of the function: for an obscure reason, MATLAB cannot 
% read the PET RTstruct of patient 40. The problem is solved by first
% computing the sData of the CT scan with contours, computing the sData of
% the PET scan without contours, and then downsampling the contours of the 
% sData of CT into the sData PET. 
% -------------------------------------------------------------------------
% INPUTS:
% - pathPET: Full path to the directory hosting the HN DICOM data of the 
%            PET scan of patient 40.
% - pathCT:  Full path to the directory hosting the HN DICOM data of the 
%            CT scan of patient 40.
% -------------------------------------------------------------------------
% OUTPUTS:
% - sDataPET: sData file of the PET scan.
% - sDataCT : sData file of the CT scan.
% -------------------------------------------------------------------------
% AUTHOR(S): Martin Vallieres <mart.vallieres@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: July 2015
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



% READING AND PROCESSING (CONTOURS) THE CT SCAN
[sDataCT] = readDICOMdir(pathCT,0);


% READING AND PROCESSING (CONTOURS) THE PET SCAN
waitB = 0;
dicomPath = pathPET;
elements = dir(dicomPath);
nElements = length(elements);
volume = cell(1,1,nElements);
dicomHeaders = [];
RTstruct = [];
REG = [];
sliceNumber = 0;
for elementNumber = 1:nElements
    elementName = elements(elementNumber).name;
    if ~strcmp(elementName,'.') && ~strcmp(elementName,'..') % Good enough for Linux, add conditions for MAC and Windows.
        elementFullFile = fullfile(dicomPath,elementName);
        if isdicom(elementFullFile)
            try
                tmp = dicominfo(elementFullFile);
            end
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
slicePositions = zeros(1,nSlices);
for sliceNumber = 1:nSlices
    slicePositions(sliceNumber) = dicomHeaders(sliceNumber).ImagePositionPatient(index);
end
[~,indices] = sort(slicePositions);
volume = cell2mat(volume(indices));
dicomHeaders = dicomHeaders(indices);
sDataPET = cell(1,5);
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
sDataPET{1} = struct('Cell_1','Explanation of cell content', ...
                  'Cell_2','Imaging data and ROI defintion (if applicable)', ...
                  'Cell_3','DICOM headers of imaging data', ...
                  'Cell_4','DICOM RTstruct (if applicable)', ...
                  'Cell_5','DICOM REGstruct (if applicable)');

sDataPET{2}.scan.volume = volume;
sDataPET{2}.scan.orientation = orientation;
try sDataPET{2}.scan.pixelW = dicomHeaders(1).PixelSpacing(1); catch sData{2}.scan.pixelW = []; end
try sDataPET{2}.scan.sliceT = dicomHeaders(1).SliceThickness; catch sData{2}.scan.sliceT = []; end
s1 = round(0.5*nSlices); s2 = round(0.5*nSlices) + 1; % Slices selected to calculate slice spacing
sDataPET{2}.scan.sliceS = sqrt(sum((dicomHeaders(s1).ImagePositionPatient - dicomHeaders(s2).ImagePositionPatient).^2)); % Slice Spacing
sDataPET{2}.type = type;
sDataPET{3} = dicomHeaders;
sDataPET{4} = RTstruct;
sDataPET{5} = REG;


% ADJUSTING PET CONTOURS ACCORDING TO CT CONTOURS.
sDataPET{2}.scan.contour = sDataCT{2}.scan.contour;
nContour = length(sDataPET{2}.scan.contour);
downF = sDataPET{2}.scan.pixelW/sDataCT{2}.scan.pixelW;
offsetStart = round(abs(sDataPET{3}(1).ImagePositionPatient(2)-sDataCT{3}(1).ImagePositionPatient(2))*sDataCT{2}.scan.pixelW);
offsetEnd = round((size(sDataCT{2}.scan.volume,1)*sDataCT{2}.scan.pixelW - size(sDataPET{2}.scan.volume,1)*sDataPET{2}.scan.pixelW)/sDataCT{2}.scan.pixelW)-offsetStart;
for i = 1:nContour
    sDataPET{2}.scan.contour(i).boxBound(1:2,1) = round((sDataCT{2}.scan.contour(i).boxBound(1:2,1) - offsetStart)/downF); sDataPET{2}.scan.contour(i).boxBound(1:2,2) = round((sDataCT{2}.scan.contour(i).boxBound(1:2,2) - offsetEnd)/downF);
    boxBound = sDataPET{2}.scan.contour(i).boxBound;
    nSlices = size(sDataPET{2}.scan.contour(i).boxMask,3);
    sDataPET{2}.scan.contour(i).boxMask = zeros(boxBound(1,2)-boxBound(1,1)+1,boxBound(2,2)-boxBound(2,1)+1,nSlices);
    for j = 1:nSlices
        sDataPET{2}.scan.contour(i).boxMask(:,:,j) = imresize(sDataCT{2}.scan.contour(i).boxMask(:,:,j),[boxBound(1,2)-boxBound(1,1)+1,boxBound(2,2)-boxBound(2,1)+1],'nearest');
    end
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