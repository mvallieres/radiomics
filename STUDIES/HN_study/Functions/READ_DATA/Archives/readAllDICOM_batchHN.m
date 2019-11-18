function readAllDICOM_batchHN(pathDICOM,nPatient,nBatch,matlabPATH)
% -------------------------------------------------------------------------
% function readAllDICOM_HN(pathDICOM,nPatient)
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
% - nPatient: Number of patients to read.
% - nBatch: Number of parallel batch.
% - matlabPATH: Full path to the MATLAB excutable on the system.
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
cd(pathDICOM), cd .., system(['mv ',pathDICOM,' DICOM']); cd('DICOM'), pathDICOM = pwd;
cd .., mkdir('DATA'), cd('DATA'), pathDATA = pwd; cd(pathDICOM)
mkdir('SCRIPTS'), cd('SCRIPTS'), pathScripts = pwd;
time = 60; % Number of seconds to wait before checking if parallel computations are done

% PRODUCE BATCH COMPUTATIONS
[patients] = batchPatients(nPatient,nBatch);
save('workspace'), pause(5);
for i = 1:nBatch
    nameScript = ['Script_readDICOM',num2str(i),'.m'];
    fid = fopen(nameScript,'w');
    fprintf(fid,'load(''workspace'')\n');
    fprintf(fid,['readAllDICOM_HN(pathDICOM,patients{',num2str(i),'})\n']);
    fprintf(fid,'clear all');
    fclose(fid);
    system([matlabPATH,' -nojvm -nodisplay -nodesktop -nosplash < ',nameScript,' >& ',nameScript(1:end-1),'log &']);
end

% WAITING LOOP
cd(pathDATA)
while 1
    pause(time);
    check = zeros(nPatient,1);
    for i = 1:nPatient
        check(i) = exist(['Patient',num2str(i),'_CT.mat']);
    end
    if sum(check) == nPatient*2
        break
    end
end

% SOLVING OTHER BUGS (the RTstruct of Patient24_PET, Patient30_CT and Patient51_CT was not read properly: No GTV included)
[sData] = correctPatient24_HN(pathDATA); cd(pathDATA), save('Patient24_PET','sData')
[sData] = correctPatient30_HN(pathDATA); cd(pathDATA), save('Patient30_CT','sData')
[sData] = correctPatient51_HN(pathDATA); cd(pathDATA), save('Patient51_CT','sData')

cd(pathScripts)
delete('workspace.mat')
cd(startpath)
end