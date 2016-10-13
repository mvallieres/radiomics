function [SUVmap] = computeSUVmap(rawPET,dicomH)
% -------------------------------------------------------------------------
% function [SUVmap] = computeSUVmap(rawPET,dicomH)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function computes the SUVmap of a raw input PET volume. It is 
% assumed that the calibration factor was applied beforehand to the PET 
% volume (e.g., rawPET = rawPET*RescaleSlope + RescaleIntercept).
% -------------------------------------------------------------------------
% INPUTS:
% - rawPET: 3D array representing the PET volume in raw format.
% - dicomH: DICOM header of one of the corresponding slice of 'rawPET'.
% -------------------------------------------------------------------------
% OUTPUTS:
% - SUVmap: 'rawPET' converted to SUVs (standard uptake values).
% -------------------------------------------------------------------------
% AUTHOR(S): 
% - Martin Vallieres <mart.vallieres@gmail.com>
% - Issam El Naqa <ielnaqa@med.umich.edu>
% - CERR development team <http://www.cerr.info/>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: May 2015
%--------------------------------------------------------------------------
% STATEMENT:
% This file is part of <https://github.com/mvallieres/radiomics/>, 
% a package providing MATLAB programming tools for radiomics analysis.
% --> Copyright (C) 2015  Martin Vallieres, Issam El Naqa
% --> Copyright 2010, Joseph O. Deasy, on behalf of the CERR development team
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

% Get patient weight
if isfield(dicomH,'PatientWeight')
    weight = dicomH.PatientWeight*1000;  % in grams
elseif isfield(dicomH,'PatientsWeight')
    weight = dicomH.PatientsWeight*1000;  % in grams
else
    weight = [];
end
if isempty(weight)
    weight = 75000; % Estimation
end

try
% Get Scan time
scantime = dcm_hhmmss(dicomH.AcquisitionTime);
% Start Time for the Radiopharmaceutical Injection
injection_time = dcm_hhmmss(dicomH.RadiopharmaceuticalInformationSequence.Item_1.RadiopharmaceuticalStartTime);
% Half Life for Radionuclide
half_life = dicomH.RadiopharmaceuticalInformationSequence.Item_1.RadionuclideHalfLife;
% Total dose injected for Radionuclide
injected_dose = dicomH.RadiopharmaceuticalInformationSequence.Item_1.RadionuclideTotalDose;

% Calculate decay
decay = exp(-log(2)*(scantime-injection_time)/half_life);
% Calculate the dose decayed during procedure
injected_dose_decay = injected_dose*decay; % in Bq

catch % Estimation
    decay = exp(-log(2)*(1.75*3600)/6588); % 90 min waiting time, 15 min preparation
    injected_dose_decay = 420000000 * decay; % 420 MBq
end

% Calculate SUV
SUVmap = rawPET*weight/injected_dose_decay;

end

% CERR UTILITY FUNCTION (can be found at: https://github.com/adityaapte/CERR)
function [totSec] = dcm_hhmmss(dateStr)
if ~ischar(dateStr)
    dateStr = num2str(dateStr);
end
hh = str2double(dateStr(1,1:2));
mm = str2double(dateStr(1,3:4));
ss = str2double(dateStr(1,5:6));
totSec = hh*60*60 + mm*60 + ss;
end