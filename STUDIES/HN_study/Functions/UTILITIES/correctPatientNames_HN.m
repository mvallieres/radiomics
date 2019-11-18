function correctPatientNames_HN(pathDATA,cohortID,patients)
% -------------------------------------------------------------------------
% function correctPatientNames_HN(pathDATA,cohortID,patients)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% Corrects the patient names to distinguish between CTsim and CT of PET-CT
% scans.
% -------------------------------------------------------------------------
% INPUTS:
% 1. pathDATA: Full path to the patient data directory
%              --> Ex: /myProject/WORKSPACE/DATA
% 2. cohortID: String specifying the cohort name (first part of patient ID)
%              --> Ex: 'HGJ'
% 3. patients: Vector of numerical value specifying patients to correct in 
%              the (cohortID) cohort.
%              --> Ex: 1:89
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

startpath = pwd;
cd(pathDATA)
nPatient = numel(patients);

for i = 1:nPatient
    if ~exist([cohortID,'_',num2str(patients(i),'%.3i'),'_CT.CTscan.mat'],'file') % If it exists, CTsim is also the CT of PET-CT, and all is ok my friend. If it does not exist, we have a 'CT1' and a 'CT2' that we need to differentiate.
        load([cohortID,'_',num2str(patients(i),'%.3i'),'_CT1.CTscan.mat']), sDataCT1 = sData;
        load([cohortID,'_',num2str(patients(i),'%.3i'),'_CT2.CTscan.mat']), sDataCT2 = sData;
        load([cohortID,'_',num2str(patients(i),'%.3i'),'_PT.PTscan.mat']), sDataPT = sData; clear sData;
        classPT = sDataPT{3}(1).FrameOfReferenceUID;
        classCT1 = sDataCT1{3}(1).FrameOfReferenceUID;
        classCT2 = sDataCT2{3}(1).FrameOfReferenceUID;
        if strcmp(classCT1,classCT2) % Need to dfferentiate by SOPClassUID
            classPT = sDataPT{3}(1).SOPClassUID;
            classCT1 = sDataCT1{3}(1).SOPClassUID;
            classCT2 = sDataCT2{3}(1).SOPClassUID;
        end
        if strcmp(classCT1,classPT)
            system(['mv ',cohortID,'_',num2str(patients(i),'%.3i'),'_CT1.CTscan.mat',' ',cohortID,'_',num2str(patients(i),'%.3i'),'_CT.CTscan.mat']);
            system(['mv ',cohortID,'_',num2str(patients(i),'%.3i'),'_CT2.CTscan.mat',' ',cohortID,'_',num2str(patients(i),'%.3i'),'_CTsim.CTscan.mat']);
        elseif strcmp(classCT2,classPT)
            system(['mv ',cohortID,'_',num2str(patients(i),'%.3i'),'_CT2.CTscan.mat',' ',cohortID,'_',num2str(patients(i),'%.3i'),'_CT.CTscan.mat']);
            system(['mv ',cohortID,'_',num2str(patients(i),'%.3i'),'_CT1.CTscan.mat',' ',cohortID,'_',num2str(patients(i),'%.3i'),'_CTsim.CTscan.mat']);
        elseif strcmp(classCT1,classCT2) % some CHUM patients needs special care
            if ~isempty(strfind(sDataCT1{3}(1).SeriesDescription,'kVCT'))
                system(['mv ',cohortID,'_',num2str(patients(i),'%.3i'),'_CT2.CTscan.mat',' ',cohortID,'_',num2str(patients(i),'%.3i'),'_CT.CTscan.mat']);
                system(['mv ',cohortID,'_',num2str(patients(i),'%.3i'),'_CT1.CTscan.mat',' ',cohortID,'_',num2str(patients(i),'%.3i'),'_CTsim.CTscan.mat']);
            elseif ~isempty(strfind(sDataCT2{3}(1).SeriesDescription,'kVCT'))
                system(['mv ',cohortID,'_',num2str(patients(i),'%.3i'),'_CT1.CTscan.mat',' ',cohortID,'_',num2str(patients(i),'%.3i'),'_CT.CTscan.mat']);
                system(['mv ',cohortID,'_',num2str(patients(i),'%.3i'),'_CT2.CTscan.mat',' ',cohortID,'_',num2str(patients(i),'%.3i'),'_CTsim.CTscan.mat']);
            end
        elseif strcmp(cohortID,'CHUM') && patients(i) == 4 % Special care for this patient
            sDataCT1([4:6]) = sDataPT([4:6]); sDataCT1{2}.scan.contour = sDataPT{2}.scan.contour;
            sData = sDataCT1; save('CHUM_004_CT.CTscan.mat','sData','-v7.3');
            sData = sDataCT2; save('CHUM_004_CTsim.CTscan.mat','sData','-v7.3');
            system('rm CHUM_004_CT1.CTscan.mat');
            system('rm CHUM_004_CT2.CTscan.mat');
        end
    end
    if strcmp(cohortID,'CHUS') && patients(i)==35 % Remove the comma in the ROI name of the GTV node --> very specific implementation that should be solved during the creation of the RTstruct
        load('CHUS_035_PT.PTscan.mat')
        [cont] = findContours_HN(sData,{',GTV2'},'PET');
        name = sData{2}.scan.contour(cont).name;
        indComma = strfind(name,','); name(indComma) = [];
        sData{2}.scan.contour(cont).name = name;
        save('CHUS_035_PT.PTscan.mat','sData','-v7.3');
        load('CHUS_035_CT.CTscan.mat')
        [cont] = findContours_HN(sData,{',GTV2'},'CT');
        name = sData{2}.scan.contour(cont).name;
        indComma = strfind(name,','); name(indComma) = [];
        sData{2}.scan.contour(cont).name = name;
        save('CHUS_035_CT.CTscan.mat','sData','-v7.3');
    end
end

cd(startpath)
end