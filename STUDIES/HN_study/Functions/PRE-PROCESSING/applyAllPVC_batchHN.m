function applyAllPVC_batchHN(pathWORK,patientNames,nBatch,matlabPATH)
% -------------------------------------------------------------------------
% function  applyAllPVC_batchHN(pathWORK,nPatient,nBatch,matlabPATH)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function applies partial-volume corrections for all PET volume of
% the head and neck cohort according to the methodology developed in ref.
% [1].
% -------------------------------------------------------------------------
% REFERENCE:
% [1] Boussion, N. et al. (2009). Incorporation of wavelet-based denoising
%     in iterative deconvolution for partial volume correction in 
%     whole-body PET imaging. Eur J Nucl Med Mol Imaging, 36(7), 1064-1075.
% -------------------------------------------------------------------------
% INPUTS:
% - pathWORK: Full path to the HN WORKSPACE directory.
% - patientNames: Cell of strings specifying the patient names to read.
% - nBatch: Number of parallel batch.
% - matlabPATH: Full path to the MATLAB excutable on the system.
% -------------------------------------------------------------------------
% AUTHOR(S): Martin Vallieres <mart.vallieres@gmail.com>
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

startpath = pwd;
cd([pathWORK,'/DATA'])
mkdir('batchLog_PVC'), cd('batchLog_PVC'), pathBatch = pwd;
time = 60; % Number of seconds to wait before checking if parallel computations are done

% PRODUCE BATCH COMPUTATIONS
[temp] = batchPatients(numel(patientNames),nBatch);
patients = cell(1,numel(temp));
for i = 1:numel(temp)
    patients{i} = cell(numel(temp{i}),1);
    patients{i}(:) = patientNames(temp{i}(:));
end
cd(pathBatch), save('workspace','pathWORK','patients','pathBatch'), pause(5);
for i = 1:nBatch
    nameScript = ['batch',num2str(i),'_script.m'];
    fid = fopen(nameScript,'w');
    fprintf(fid,'dummy = 0;\n');
    fprintf(fid,'load(''workspace'')\n');
    fprintf(fid,['applyAllPVC_HN(pathWORK,patients{',num2str(i),'})\n']);
    fprintf(fid,'cd(pathBatch)\n');
    fprintf(fid,['save(''batch',num2str(i),'_end'',''dummy'')\n']);
    fprintf(fid,'clear all');
    fclose(fid);
    system([matlabPATH,' -nojvm -nodisplay -nodesktop -nosplash < ',nameScript,' >& ',nameScript(1:end-1),'log &']);
end

% WAITING LOOP
waitBatch(pathBatch,time,nBatch)
cd(pathBatch)
delete('workspace.mat')

cd(startpath)
end