function organizeRF_LGG(pathTCGA,pathSTUDY,pathTextModels,pathVasModels,pathRF)
% -------------------------------------------------------------------------
% function organizeRF_LGG(pathTCGA,pathSTUDY,pathTextModels,pathVasModels,pathRF)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function organizes the data for the random forest experiments of this study.
% -------------------------------------------------------------------------
% OUTPUTS: MATLAB file named 'dataStruct.mat' in /myProject/WORKSPACE/RANDOM_FORESTS
% -------------------------------------------------------------------------
% AUTHOR(S): 
% - Martin Vallieres <mart.vallieres@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: January 2017
%--------------------------------------------------------------------------
% STATEMENT:
% This file is part of <https://github.com/mvallieres/radiomics/>, 
% a package providing MATLAB programming tools for radiomics analysis.
% --> Copyright (C) 2015-2017  Martin Vallieres
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

cd(pathSTUDY), load('outcomes'), load('lesionNum')
cd(pathTCGA), load('clinical')
skipFields = {{'IDHsub','primaryRT','primaryPharma','sample','RseqCluster','MethylationCluster','miRCluster','CNCluster','RPPACluster','OncosignCluster','COCCluster','ATRX','CIC','EGFR','TP53'},{'IDHsub','primaryRT','primaryPharma','sample','RseqCluster','MethylationCluster','miRCluster','CNCluster','RPPACluster','OncosignCluster','COCCluster','ATRX','CIC','EGFR','TP53'},{'sample','RseqCluster','MethylationCluster','miRCluster','CNCluster','RPPACluster','OncosignCluster','COCCluster','ATRX','CIC','EGFR','TP53'},{'grade','primaryRT','primaryPharma','sample','RseqCluster','MethylationCluster','miRCluster','CNCluster','RPPACluster','OncosignCluster','COCCluster','ATRX','CIC','EGFR','TP53'}}; % One cell for each outcome (IDH1,IDHcodel,progression,lowGrade)
vasariSets = {'VASARI','VASARI','VASARI','VASARI'}; % One for each outcome (IDH1,IDHcodel,progression,lowGrade)
textureSets = {'T1CE_T2W','T1CE_T2W','T1W_T2W','T1CE_T2W'}; % One for each outcome (IDH1,IDHcodel,progression,lowGrade)

nameOutcomes = fieldnames(outcomes); nOutcomes = numel(nameOutcomes);
nameClinical = fieldnames(clinical.data); nClinical = numel(nameClinical);

dataStruct = struct;
for o = 1:nOutcomes
    nameOutcome = nameOutcomes{o};
    dataStruct.(nameOutcome) = struct;
    lesions = lesionNum.(nameOutcome);
    outcome = outcomes.(nameOutcome);
    
    % GO GET VASARI
    cd(fullfile(pathVasModels,nameOutcome,vasariSets{o})), load('finalModel')
    vasari = finalModel.Data; clear finalModel
    nVas = size(vasari,2);
    
    % GO GET TEXT
    cd(fullfile(pathTextModels,nameOutcome,textureSets{o})), load('finalModel')
    text = finalModel.Data; clear finalModel
    nText = size(text,2);
    
    % PROCESS VASARI AND TEXT
    strString = []; cat = [];
    for i = 1:nVas
        feat = ['V',num2str(i)];
        eval([feat,' = vasari(:,i);']);
        strString = [strString,feat,','];
        cat = [cat,0];
    end
    for i = 1:nText
        feat = ['T',num2str(i)];
        eval([feat,' = text(:,i);']);
        strString = [strString,feat,','];
        cat = [cat,0];
    end
    
    % PROCESSING CLINICAL
    for i = 1:nClinical
        feat = nameClinical{i};
        if ~sum(strcmp(feat,skipFields{o}))
            eval([feat,' = clinical.data.',feat,'(lesions);']);
            strString = [strString,feat,','];
            if iscell(clinical.data.(feat))
                cat = [cat,1];
            else
                cat = [cat,0];
            end
        end
    end
    strString(end) = [];
    
    % ARRANGE DATA STRUCT
    eval(['tab = table(',strString,');']);
    dataStruct.(nameOutcome).data = tab;
    dataStruct.(nameOutcome).categories = cat;
    dataStruct.(nameOutcome).outcome = outcome;
    dataStruct.(nameOutcome).nVas = nVas;
    dataStruct.(nameOutcome).nText = nText;
end

cd(pathRF), save('dataStruct','dataStruct') % Final structure organizing all data for random forests experiments.

cd(startpath)
end