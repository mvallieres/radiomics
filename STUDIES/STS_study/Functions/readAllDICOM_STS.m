function readAllDICOM_STS(pathDICOM,nPatient)
% -------------------------------------------------------------------------
% function readAllDICOM_STS(pathDICOM,nPatient)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% Reads all soft-tissue sarcoma (STS) DICOM imaging data downloaded from 
% The Cancer Imaging Archive (TCIA) website at: 
% <http://dx.doi.org/10.7937/K9/TCIA.2015.7GO2GSKS>. Imaging data is then 
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
% [1] Vallieres, M. et al. (2015). A radiomics model from joint FDG-PET and 
%     MRI texture features for the prediction of lung metastases in soft-tissue 
%     sarcomas of the extremities. Physics in Medicine and Biology, 60(14), 
%     5471-5496. doi:10.1088/0031-9155/60/14/5471
% -------------------------------------------------------------------------
% INPUTS:
% - pathDICOM: Full path to the directory hosting the STS DICOM data 
%              downloaded from the TCIA website.
% - nPatient: Number of patients to read. 
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


waitbar = 0;
cd(pathDICOM), cd .., mkdir DATA, cd('DATA'); pathDATA=pwd; cd(pathDICOM);

% ORGANIZING TCIA DATA
cd .., system(['mv ',pathDICOM,' DICOM']); cd('DICOM'); pathDICOM=pwd;
names = {'CT','PET','T1','T2FS','T1reg','T2FSreg'};
nFolder = numel(names);
fprintf('\n')
for i = 1:nPatient
    fprintf('ORGANIZING TCIA DATA FOR PATIENT %u ... ',i)
    cd(['STS_',num2str(i,'%0.3i')]), pathPatient = pwd;
    level1 = dir;
    for f = 1:nFolder
        mkdir(names{f});
    end
    for j = 1:numel(level1) - 2
        cd(level1(j+2).name), level2 = dir;
        for k = 1:numel(level2) - 2
            cd(level2(k+2).name)
            info = dicominfo('000000.dcm');
            if ~isempty(strfind(info.SeriesDescription,'PET')) && (~isempty(strfind(info.SeriesDescription,'T1')) || ~isempty(strfind(info.SeriesDescription,'T!')))
                folder = 5;
            elseif ~isempty(strfind(info.SeriesDescription,'PET')) && (~isempty(strfind(info.SeriesDescription,'T2FS')) || ~isempty(strfind(info.SeriesDescription,'STIR')))
                folder = 6;
            elseif ~isempty(strfind(info.StudyDescription,'PET')) && isempty(strfind(info.SeriesDescription,'T1')) && isempty(strfind(info.SeriesDescription,'T2FS')) && isempty(strfind(info.SeriesDescription,'STIR'))
                if strcmp(info.Modality,'RTSTRUCT')
                    if ~isempty(strfind(info.SeriesDescription,'PET'))
                        folder = 2;
                    elseif ~isempty(strfind(info.SeriesDescription,'CT'))
                        folder = 1;
                    end
                else
                    if strcmp(info.Modality,'PT')
                        folder = 2;
                    elseif strcmp(info.Modality,'CT')
                        folder = 1;
                    end
                end
            elseif ~isempty(strfind(info.SeriesDescription,'T1')) || ~isempty(strfind(info.SeriesDescription,'t1'))
                folder = 3;
            else
                folder = 4;
            end
            content = dir;
            for c =1:numel(content)
                if strfind(content(c).name,'.dcm')
                    if strcmp(info.Modality,'RTSTRUCT')
                        nameSave = 'RTstruct.dcm';
                    else
                        nameSave = content(c).name;
                    end
                    system(['mv ',content(c).name,' ',pathPatient,'/',names{folder},'/',nameSave]);
                end
            end
            cd ..
        end
        cd .., system(['rm -r ',level1(j+2).name]);
    end
    cd ..
    fprintf('DONE\n')
end
cd(pathDICOM)
fid = fopen('README.txt','w');
fprintf(fid,'T1reg: T1 scan registered and resampled (axially) onto PET\n');
fprintf(fid,'T2FSreg: T2FS scan registered and resampled (axially) onto PET');
fclose(fid);


% READING AND SAVING ALL DATA
fprintf('\n')
names = {'PET','T1','T2FS','T1reg','T2FSreg'};
nFolder = numel(names);
for i = 1:nPatient
    namePatient = ['Patient',num2str(i)];
    for j = 1:nFolder
        cd([pathDICOM,'/','STS_',num2str(i,'%0.3i'),'/',names{j}])
        fprintf(['PROCESSING ',namePatient,'_',names{j},' ... '])
        [sData] = readDICOMdir(pwd,waitbar);
        cd(pathDATA), save([namePatient,'_',names{j}],'sData')
        fprintf('DONE\n')
    end
    fprintf('\n')
end


% SOLVING KNOWN BUGS (1. Slice thickness of MRI scans of Patient 47 not found in DICOM headers coming out of MIM)
cd(pathDATA)
load('Patient47_T1'), sData{2}.scan.sliceT = 4; save('Patient47_T1','sData')
load('Patient47_T2FS'), sData{2}.scan.sliceT = 4; save('Patient47_T2FS','sData')

cd ..
end