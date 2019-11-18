function readAllDICOM_HN(pathDICOM,patients)
% -------------------------------------------------------------------------
% function readAllDICOM_HN(pathDICOM,patients)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% Reads all head and neck (HN) DICOM imaging data downloaded from 
% The Cancer Imaging Archive (TCIA) website at: 
% <http://dx.doi.org/xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx>. Imaging data is then 
% organized in 'sData' format under a folder 'DATA'. 'sData' files are 
% cells of structures organizing the content of the volume data, DICOM 
% headers, DICOM RTstruct* (used to compute the ROI) and  DICOM REGstruct* 
% (used to register a MRI volume to a PET volume).
%   * If present in the directory
%   --> sData{1}: Explanation of cell content
%   --> sData{2}: Imaging data and ROI defintion (if applicable)
%   --> sData{3}: DICOM headers of imaging data
%   --> sData{4}: DICOM RTstruct (if applicable)
%   --> sData{5}: DICOM REGstruct (if applicable)
% -------------------------------------------------------------------------
% REFERENCE:
% [1] Vallieres, M. et al. (2015). FDG-PET/CT radiomics models for the 
%     early prediction of different tumour outcomes in head and neck cancer.
%     The Journal of Nuclear Medicine, aa(bb), xxx-yyy. 
%     doi:
% -------------------------------------------------------------------------
% INPUTS:
% - pathDICOM: Full path to the directory hosting the HN DICOM data 
%              downloaded from the TCIA website.
% - patients: Vector of patient number to analyze.
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

startpath = pwd;
waitbar = 0;
cd(pathDICOM), cd .., cd('DATA'); pathDATA = pwd; cd(pathDICOM);
nPatient = length(patients);

% ORGANIZING TCIA DATA
% Put code here once data is online on TCIA website%
% Temporary code
for i = 1:nPatient
    try
        cd(['Patient',num2str(patients(i))])
        cd ..
        system(['mv Patient',num2str(patients(i)),' HN_',num2str(patients(i),'%0.3i')]);
    end
end


% READING AND SAVING ALL DATA
fprintf('\n')
names = {'PET','CT'};
nFolder = numel(names);
for i = 1:nPatient
    namePatient = ['Patient',num2str(patients(i))];
    if patients(i) == 40 % solving known bug: MATLAB cannot read the PET RTstruct of patient 40
        fprintf(['PROCESSING ',namePatient,' (known bug) ... '])
        pathPET = [pathDICOM,'/','HN_',num2str(patients(i),'%0.3i'),'/PET'];
        pathCT = [pathDICOM,'/','HN_',num2str(patients(i),'%0.3i'),'/CT'];
        [sDataPET,sDataCT] = readPatient40_HN(pathPET,pathCT);
        sData = sDataPET; cd(pathDATA), save([namePatient,'_PET'],'sData')
        sData = sDataCT; cd(pathDATA), save([namePatient,'_CT'],'sData')
        fprintf('DONE\n')
    else
        for j = 1:nFolder
            cd([pathDICOM,'/','HN_',num2str(patients(i),'%0.3i'),'/',names{j}])
            fprintf(['PROCESSING ',namePatient,'_',names{j},' ... '])
            [sData] = readDICOMdir(pwd,waitbar);
            cd(pathDATA), save([namePatient,'_',names{j}],'sData')
            fprintf('DONE\n')
        end
    end
    fprintf('\n')
end

cd(startpath)
end