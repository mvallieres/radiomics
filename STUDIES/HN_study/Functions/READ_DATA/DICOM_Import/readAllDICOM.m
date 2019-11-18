function readAllDICOM(pathRead,pathSave,nBatch,nameSaveOption)
% -------------------------------------------------------------------------
% function readAllDICOM(pathRead,pathSave,batchFlag,nameSaveOption)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function reads the DICOM content of all the sub-folder tree of a 
% starting directory defined by 'pathRead'. It then organizes the data 
% (files throughout the starting directory are associated by 
% 'SeriesInstanceUID') in a cell of structures called 'sData' (stands for 
% STAMP data), and it finally computes the region of interest (ROI) defined
% by an associated RTstruct, and save as well all the  REG, RTdose and 
% RTplan structures(if present in the sub-folder tree of the starting 
% directory). 
% All sData structures hereby created are saved in 'pathSave' with a name 
% varying with the variable 'nameSaveOption'. 
%
% DIFFERENTIATION/ASSOCIATION OF DICOM FILES: 1)imaging, 2)RTstruct, 3)REG, 4)RTdose, 5)RTplan. 
% 1) Imaging volumes are differentiated by the 'SeriesInstanceUID' field
% 2) Association between a RTstruct and an imaging volume is performed with:
%    - 'ReferencedFrameOfReferenceSequence.Item_1.FrameOfReferenceUID' field, or
%    - 'FrameOfReferenceUID' field. 
% 3) Association between a REG and an imaging volume is performed with:
%    - 'RegistrationSequence.Item_2.FrameOfReferenceUID' field, or
%    - 'FrameOfReferenceUID' field. 
% 4) Association between a RTdose and an imaging volume is performed with:
%    - 'FrameOfReferenceUID' field.
% 5) Association between a RTplan and an imaging volume is performed with:
%    - 'FrameOfReferenceUID' field.
%
% If multiple imaging volumes have the same 'FrameOfReferenceUID' field, 
% they will all be assigned a RTstruct, REG, RTdose and RTplan possessing 
% that same field.
% 
% IMPORTANT: This function has only been tested on Ubuntu 14.04. May work 
%            using MAC, will fail using Windows.
% -------------------------------------------------------------------------
% INPUTS:
% 1. pathRead: String specifying the full path to the starting directory 
%              where the DICOM path to read are located.
%              --> Ex: '/home/myStudy/DICOM'
% 2. pathSave: String specifying the full path to the directory where to 
%              save all the sData structures created by the current function.
%              --> Ex: '/home/myStudy/DATA'
% 3. nBatch: Numerical value specifying the number of batch to use in the
%            parallel computations (use 0 for serial).
%            --> Ex: 8
% 4. nameSaveOption: (optional). If this argument is not present (default)
%                    all sData structures will be saved with a filename
%                    defined by the 'Modality' field present in the
%                    DICOM header of their corresponding imaging volume. If
%                    multiple volumes of the same modality are present for 
%                    the same 'PatientID', the different volumes will be 
%                    numerated (e.g. CT1, CT2, CT3, etc.). If this argument 
%                    is present and set to 'folder', all sData structures 
%                    will be saved with a filename defined by the name of 
%                    the folder in which all DICOM files of a given imaging 
%                    volume are located (useful to keep track of meaningful 
%                    names). The 'folder' option assumes that all imaging 
%                    DICOM files of a given imaging volume are present in 
%                    the same directory. Note: the 'folder' argument may 
%                    lead to overwritting of data if not organized properly.
%                    Finally, if this argument is present and set to 
%                    'modality',  all sData structures will be saved with a 
%                    filename defined by the 'Modality' field in the DICOM
%                    header, but without enumerating multiple volumes of
%                    the same modality and same 'PatientID' (this may lead
%                    to overwritting, but is faster if the user is sure
%                    that only one scan of each modality for each patient
%                    is present in 'pathRead'.
%                    --> Options: - No argument
%                                 - 'folder'
%                                 - 'modality'  
% -------------------------------------------------------------------------
% OUTPUTS: Multiple sData files
% - sData: Cell of structures organizing the content of the volume data, 
%          DICOM headers, DICOM RTstruct* (used to compute the ROI) and 
%          DICOM REGstruct* (used to register a given volume to another),
%          DICOM RTdose* and DICOM RTplan*.
%          * If present in the sub-folder tree
%    --> sData{1}: Explanation of cell content
%    --> sData{2}: Imaging data and ROI defintion
%    --> sData{3}: DICOM headers of imaging data
%    --> sData{4}: DICOM RTstruct (if applicable)
%    --> sData{5}: DICOM REGstruct (if applicable)
%    --> sData{6}: DICOM RTdose (if applicable)
%    --> sData{7}: DICOM RTplan (if applicable)
% -------------------------------------------------------------------------
% AUTHOR(S): 
% - Martin Vallieres <mart.vallieres@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: February 2016
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
warning off

% INITIALIZATION
cd(pathRead)
stackPathRS = {}; stackFrameRS = {}; % Full path and FrameOfReferenceUID of the RTstruct files.
stackPathREG = {}; stackFrameREG = {}; % Full path and FrameOfReferenceUID of the REG files.
stackPathRD = {}; stackFrameRD = {}; % Full path and FrameOfReferenceUID of the RTdosefiles.
stackPathRP = {}; stackFrameRP = {}; % Full path and FrameOfReferenceUID of the RTplan files.
cellSeriesID = {}; % Cell of 'SeriesInstanceUID' of different imaging volumes (string).
cellFrameID = {}; % Cell of 'FrameOfReferenceUID' of different imaging volumes (string). Cell index is associated to the index of cellSeriesID.
cellPathImages = {}; % Cell of paths to the different imaging volumes. Cell index is associated to the index of cellSeriesID. Each cell in turn contains all the different path to the dicom images of a given volume.
cellPathRS = {}; % Cell of paths to the different RTstruct file (struct). Cell index is associated to the index of cellSeriesID.
cellPathREG = {}; % Cell of paths to the different REG file (struct). Cell index is associated to the index of cellSeriesID.
cellPathRD = {}; % Cell of paths to the different RTdose file (struct). Cell index is associated to the index of cellSeriesID.
cellPathRP = {}; % Cell of paths to the different RTplan file (struct). Cell index is associated to the index of cellSeriesID.
nameSave = {}; % Cell of filenames for all created sData structures (string). Cell index is associated to the index of cellSeriesID.
fidErr = fopen('readAllDICOM_err','w');
fidRead = fopen('readAllDICOM_log','w');

% SCANNING ALL FOLDERS IN INITIAL DIRECTORY
stackFolder = {pathRead}; % Stack of folder paths to scan
ind = strfind(stackFolder{1},'/'); nameDir = stackFolder{1}((ind(end)+1):end);
indUnder = strfind(nameDir,'_');
if ~isempty(indUnder)
    for under = numel(indUnder)
        nameDir = [nameDir(1:indUnder(under)-1),'\',nameDir(indUnder(under):end)];
    end
end
while ~isempty(stackFolder)
    cd(stackFolder{1})
    fprintf(fidRead,['\nScanning files in directory: ',stackFolder{1},' ... ']);
    list = dir; nElem = numel(list);
    for i = 1:nElem
        if ~strcmp(list(i).name,'.') && ~strcmp(list(i).name,'..') && isempty(strfind(list(i).name,'DS_Store')) % Good enough for linux, maybe add other conditions for MAC and Windows
            if list(i).isdir
                cd(list(i).name), stackFolder = [stackFolder;pwd]; cd ..
            elseif isdicom(list(i).name)
                try
                    info = dicominfo(list(i).name);
                    if strcmp(info.Modality,'MR') || strcmp(info.Modality,'PT') || strcmp(info.Modality,'CT')
                        [indSeriesID] = findUIDcellIndex(info.SeriesInstanceUID,cellSeriesID);
                        if indSeriesID > numel(cellSeriesID) % New volume
                            cellSeriesID = [cellSeriesID cell(1)]; cellSeriesID{indSeriesID} = info.SeriesInstanceUID; 
                            cellFrameID = [cellFrameID cell(1)]; cellFrameID{indSeriesID} = info.FrameOfReferenceUID;
                            cellPathImages = [cellPathImages cell(1)]; cellPathImages{indSeriesID} = {};
                            cellPathRS = [cellPathRS cell(1)]; cellPathRS{indSeriesID} = {};
                            cellPathREG = [cellPathREG cell(1)]; cellPathREG{indSeriesID} = {};
                            cellPathRD = [cellPathRD cell(1)]; cellPathRD{indSeriesID} = {};
                            cellPathRP = [cellPathRP cell(1)]; cellPathRP{indSeriesID} = {};
                            nameSave = [nameSave cell(1)];
                            if nargin > 3 && strcmp(nameSaveOption,'folder')
                                path = pwd;
                                ind = strfind(path,'/'); % NOT GOOD FOR WINDOWS
                                nameSave{indSeriesID} = path((ind(end)+1):end);
                            elseif nargin > 3 && strcmp(nameSaveOption,'modality')
                                nameSave{indSeriesID} = info.Modality;
                            else
                                try
                                    nameSave{indSeriesID} = info.SeriesDescription;
                                catch
                                    nameSave{indSeriesID} = info.SeriesInstanceUID;
                                end
                            end
                        end
                        cellPathImages{indSeriesID} = [cellPathImages{indSeriesID};fullfile(pwd,list(i).name)];
                    elseif strcmp(info.Modality,'RTSTRUCT')
                        stackPathRS = [stackPathRS;fullfile(pwd,list(i).name)];
                        try
                            frameUID = info.ReferencedFrameOfReferenceSequence.Item_1.FrameOfReferenceUID;
                        catch
                            frameUID = info.FrameOfReferenceUID;
                        end
                        stackFrameRS = [stackFrameRS;frameUID];
                    elseif strcmp(info.Modality,'REG')
                        stackPathREG = [stackPathREG;fullfile(pwd,list(i).name)];
                        try
                            frameUID = info.RegistrationSequence.Item_2.FrameOfReferenceUID;
                        catch
                            frameUID = info.FrameOfReferenceUID;
                        end
                        stackFrameREG = [stackFrameREG;frameUID];
                    elseif strcmp(info.Modality,'RTDOSE')
                        stackPathRD = [stackPathRD;fullfile(pwd,list(i).name)];
                        frameUID = info.FrameOfReferenceUID;
                        stackFrameRD = [stackFrameRD;frameUID];
                    elseif strcmp(info.Modality,'RTPLAN')
                        stackPathRP = [stackPathRP;fullfile(pwd,list(i).name)];
                        frameUID = info.FrameOfReferenceUID;
                        stackFrameRP = [stackFrameRP;frameUID];
                    end
                catch
                    fprintf(fidErr,'Error while reading: %s\n',fullfile(pwd,list(i).name));
                end
            end
        end
    end
    stackFolder(1) = [];
    fprintf(fidRead,'DONE!');
end
fclose(fidErr);
fprintf(fidRead,'\n');

% ASSOCIATING ALL RTSTRUCT TO IMAGING VOLUMES
nRS = numel(stackPathRS);
if nRS
    fprintf(fidRead,'\nASSOCIATING ALL RTSTRUCT TO IMAGING VOLUMES ... ');
    for i = 1:nRS
        [indSeriesID] = findUIDcellIndex(stackFrameRS{i},cellFrameID);
        for n = 1:numel(indSeriesID)
            cellPathRS{indSeriesID(n)} = [cellPathRS{indSeriesID(n)};stackPathRS{i}];
        end
    end
    fprintf(fidRead,'DONE!');
end

% ASSOCIATING ALL REG TO IMAGING VOLUMES
nREG = numel(stackPathREG);
if nREG
    fprintf(fidRead,'\nASSOCIATING ALL REG TO IMAGING VOLUMES ... ');
    for i = 1:nREG
        [indSeriesID] = findUIDcellIndex(stackFrameREG{i},cellFrameID);
        for n = 1:numel(indSeriesID)
            cellPathREG{indSeriesID(n)} = [cellPathREG{indSeriesID(n)};stackPathREG{i}];
        end
    end
    fprintf(fidRead,'DONE!');
end


% ASSOCIATING ALL RTdose TO IMAGING VOLUMES
nRD = numel(stackPathRD);
if nRD
    fprintf(fidRead,'\nASSOCIATING ALL RTdose TO IMAGING VOLUMES ... ');
    for i = 1:nRD
        [indSeriesID] = findUIDcellIndex(stackFrameRD{i},cellFrameID);
        for n = 1:numel(indSeriesID)
            cellPathRD{indSeriesID(n)} = [cellPathRD{indSeriesID(n)};stackPathRD{i}];
        end
    end
    fprintf(fidRead,'DONE!');
end

% ASSOCIATING ALL RTplan TO IMAGING VOLUMES
nRP = numel(stackPathRP);
if nRP
    fprintf(fidRead,'\nASSOCIATING ALL RTplan TO IMAGING VOLUMES ... ');
    for i = 1:nRP
        [indSeriesID] = findUIDcellIndex(stackFrameRP{i},cellFrameID);
        for n = 1:numel(indSeriesID)
            cellPathRP{indSeriesID(n)} = [cellPathRP{indSeriesID(n)};stackPathRP{i}];
        end
    end
    fprintf(fidRead,'DONE!');
end


% READING ALL IMAGES TO CREATE ALL sData STRUCTURES
nScans = numel(cellPathImages);
if nBatch
    if nScans < nBatch
        nBatch = nScans;
    end
    fprintf(fidRead,['\n\nREADING AND PROCESSING ALL DICOM DATA USING ',num2str(nBatch),' PARALLEL BATCH ... ']);
    cd(pathRead), mkdir('BATCH_LOG'), cd('BATCH_LOG'), pathBatch = pwd;
    [scans] = batchPatients(nScans,nBatch); % Scans for each batch
    save('workspace','scans','pathSave','cellPathImages','cellPathRS','cellPathREG','cellPathRD','cellPathRP','nameSave'), pause(2)
    for i = 1:nBatch
        nScan = numel(scans{i});
        nameScript = ['batch',num2str(i),'_script.m'];
        fid = fopen(nameScript,'w');
        fprintf(fid,'dummy = 0;\n');
        fprintf(fid,'load(''workspace'')\n');
        for j = 1:nScan
            fprintf(fid,['sDataCreation_FromDICOMpaths(pathSave,cellPathImages{scans{',num2str(i),'}(',num2str(j),')},cellPathRS{scans{',num2str(i),'}(',num2str(j),')},cellPathREG{scans{',num2str(i),'}(',num2str(j),')},cellPathRD{scans{',num2str(i),'}(',num2str(j),')},cellPathRP{scans{',num2str(i),'}(',num2str(j),')},nameSave{scans{',num2str(i),'}(',num2str(j),')})\n']);
        end
        fprintf(fid,['system(''touch batch',num2str(i),'_end'');\n']);
        fprintf(fid,'clear all');
        fclose(fid);
        system(['matlab -nojvm -nodisplay -nodesktop -nosplash < ',nameScript,' >& ',nameScript(1:end-1),'log &']);
    end
    waitBatch(pathBatch,30,nBatch)
    cd ..
    %system('rm -r BATCH_LOG');
else
    fprintf(fidRead,'\n\nREADING AND PROCESSING ALL DICOM DATA ... ');
    for i = 1:nScans
        sDataCreation_FromDICOMpaths(pathSave,cellPathImages{i},cellPathRS{i},cellPathREG{i},cellPathRD{i},cellPathRP{i},nameSave{i});
    end
end
fprintf(fidRead,'DONE!\n');


% POST-PROCESSING OF SAVED FILES IF 4th ARGUMENT IS SET TO 'modality'
if nargin < 4
    fprintf(fidRead,'\nCORRECTING NAME OF SAVED FILES ... ');
    type = {'.CTscan','.PTscan','.MRscan'}; nType = numel(type);
    cd(pathSave)
    list = dir; nElem = numel(list);
    cellNames = {};
    for i = 1:nElem
        if ~list(i).isdir
            cellNames = [cellNames,list(i).name];
        end
    end
    nNames = numel(cellNames);
    temp = '';
    for i = 1:nNames
        try
            load(cellNames{i})
            idString = sData{3}(1).PatientID;
            if strcmp(idString,temp)
                ok = 0;
            else
                ok = 1;
            end
        catch
            ok = 0;
        end
        if ok
            cellTemp = strfind(cellNames,idString);
            for c = 1:numel(cellTemp)
                if isempty(cellTemp{c})
                    cellTemp{c} = [0];
                end
            end
            idIndex = find(cell2mat(cellTemp)); 
            nIndex = numel(idIndex); nameComp = cell(1,nIndex);
            for j = 1:nIndex
                nameComp{j} = cellNames{idIndex(j)};
            end
            for j = 1:nType
                cellTemp = strfind(nameComp,type{j});
                for c = 1:numel(cellTemp)
                    if isempty(cellTemp{c})
                        cellTemp{c} = [0];
                    end
                end
                scanIndex = find(cell2mat(cellTemp)); 
                nScanIndex = numel(scanIndex);
                if nScanIndex > 0
                    for k = 1:nScanIndex
                        ind1 = numel(idString) + 1; % Getting the '_' position after the PatientID 
                        ind2 = numel(nameComp{scanIndex(k)}) - 10; % Getting the '.' position
                        name = nameComp{scanIndex(k)};
                        indWrong = strfind(name,'\');
                        if ~isempty(indWrong)
                            nInd = length(indWrong);
                            for n = 1:nInd
                                name = [name(1:indWrong(n)+n-2),'\',name(indWrong(n)+n-1:end)];
                            end
                        end
                        indWrong = strfind(name,'(');
                        if ~isempty(indWrong)
                            nInd = length(indWrong);
                            for n = 1:nInd
                                name = [name(1:indWrong(n)+n-2),'\',name(indWrong(n)+n-1:end)];
                            end
                        end
                        indWrong = strfind(name,')');
                        if ~isempty(indWrong)
                            nInd = length(indWrong);
                            for n = 1:nInd
                                name = [name(1:indWrong(n)+n-2),'\',name(indWrong(n)+n-1:end)];
                            end
                        end
                        indWrong = strfind(name,'&');
                        if ~isempty(indWrong)
                            nInd = length(indWrong);
                            for n = 1:nInd
                                name = [name(1:indWrong(n)+n-2),'\',name(indWrong(n)+n-1:end)];
                            end
                        end
                        indWrong = strfind(name,' ');
                        if ~isempty(indWrong)
                            nInd = length(indWrong);
                            for n = 1:nInd
                                name = [name(1:indWrong(n)+n-2),'\',name(indWrong(n)+n-1:end)];
                            end
                        end
                        if nScanIndex == 1
                            system(['mv ',name,' ',nameComp{scanIndex(k)}(1:ind1(1)),type{j}(2:3),nameComp{scanIndex(k)}(ind2(1):end)]);
                        else
                            system(['mv ',name,' ',nameComp{scanIndex(k)}(1:ind1(1)),type{j}(2:3),num2str(k),nameComp{scanIndex(k)}(ind2(1):end)]);
                        end
                    end
                end
            end
            temp = idString;
        end
    end
    fprintf(fidRead,'DONE!\n');
end
fclose(fidRead);

cd(startpath)
end