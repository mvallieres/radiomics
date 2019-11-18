function applyAllPVC_HN(pathWORK,patientNames)
% -------------------------------------------------------------------------
% function applyAllPVC_HN(pathWORK,patientNames)
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
nPatient = length(patientNames);

% PVE CORRECTIONS
for i = 1:nPatient
    fprintf('\n***** CORRECTING %s *****\n',patientNames{i})
    load(patientNames{i}) % Variable 'sData' now in MATLAB Workspace
    volumePVC = PVEcorrect(sData{2}.scan.volume.data,2); % 2 iterations have proved to be sufficient
    sData{2}.scan.volume.data = volumePVC;
    sData{2}.scan.info = 'PVE corrected';
    save(patientNames{i},'sData')
end

cd(startpath)
end