function organizeRadiomicsExp_LGG(pathSTUDY,pathLR,nonTextName,fSetNames,paramSEP,baselineSEP,freedomSEP,textType,textName)
% -------------------------------------------------------------------------
% function organizeRadiomicsExp_LGG(pathSTUDY,pathLR,pathCR,outcomes,timeToEvent,nonTextName,fSetNames,paramSEP,baselineSEP,freedomSEP,textType,textName)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function organizes the radiomics experiments for that study.
% -------------------------------------------------------------------------
% INPUTS:
% 1. pathSTUDY: Full path to where the organized data will be saved.
%               Defined in masterScript_LGG.m
%               --> Ex: '/myProject/WORKSPACE/STUDY_DATA'
% 2. pathLR: Full path to where the logistic regression experiments results
%            will be saved. Defined in masterScript_LGG.m
%            --> Ex: '/myProject/WORKSPACE/LOGISTIC_REGRESSION'
% 3. nonTextName: Cell of strings specifying the non-texture names used in 
%                 this study. 
%                 --> Defined in masterScript_LGG.m
% 4. fSetNames: Cell of strings specifying the names of the different
%               feature set studied.
%               --> Ex: {'T1W_T2W','T1W_T2F','T1CE_T2W','T1CE_T2F'}
% 5. paramSEP: Cell specifying all the different texture extraction
%              parameters used in this study.
%              --> Defined in masterScript_LGG.m
% 6. baselineSEP: Cell specifying all the baseline texture extraction
%                 parameters used in this study.
%                 --> Defined in masterScript_LGG.m
% 7. freedomSEP: Cell specifying all the degree of freedom on texture 
%                extraction parameters used in this study.
%                --> Defined in masterScript_LGG.m
% 8. textType: Cell of strings specifying the texture types used in this study.
%              --> Defined in masterScript_LGG.m
% 9. textName: Cell of strings specifying the texture names used in this study.s
%              --> Defined in masterScript_LGG.m
% -------------------------------------------------------------------------
% OUTPUTS: Structure named 'training' organizing all parameters for all 
%          experiments. 
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
cd(pathSTUDY), load('outcomes')

% INITIALIZATION
nNonText = numel(nonTextName);
nFset = numel(fSetNames);
nameOutcomes = fieldnames(outcomes); nOutcomes = numel(nameOutcomes);
nTextType = numel(textType);


% DATA ORGANIZATION
training = struct;
for o = 1:nOutcomes
    
    % Organizing non-texture features
    cd(pathSTUDY)
    train1 = load(['nonText_',nameOutcomes{o}]); train1 = struct2cell(train1); train1 = train1{1};
    for i = 1:nNonText
        training.(nameOutcomes{o}).nonText.(nonTextName{i}).Data = train1.(nonTextName{i}).Data;
    end

    % Organizing texture features
    for f = 1:nFset
        ind = strfind(fSetNames{f},'_');
        set1 = fSetNames{f}(1:(ind-1)); set2 = fSetNames{f}((ind+1):end);
        train1_1 = load(['text_',nameOutcomes{o},'_',set1]); train1_1 = struct2cell(train1_1); train1_1 = train1_1{1};
        train1_2 = load(['text_',nameOutcomes{o},'_',set2]); train1_2 = struct2cell(train1_2); train1_2 = train1_2{1};
        train1 = {train1_1,train1_2};
        training.(nameOutcomes{o}).text.(fSetNames{f}).param = paramSEP; 
        training.(nameOutcomes{o}).text.(fSetNames{f}).baseline = baselineSEP;
        training.(nameOutcomes{o}).text.(fSetNames{f}).freedom = freedomSEP;
        training.(nameOutcomes{o}).text.(fSetNames{f}).paramName = {'Scale','Quant.algo','Ng'};
        textCellName = {set1,set2};

        for i = 1:numel(train1)
            sizeParam = size(train1{1});
            training.(nameOutcomes{o}).text.(fSetNames{f}).(textCellName{i}) = cell(sizeParam);
            for n = 1:numel(train1{1})
                for type = 1:nTextType
                    for text = 1:numel(textName{type})
                        training.(nameOutcomes{o}).text.(fSetNames{f}).(textCellName{i}){n}.(textType{type}).(textName{type}{text}).Data = train1{i}{n}.(textType{type}).(textName{type}{text}).Data;
                    end
                end
            end
        end
    end
end

% Organizing outcomes
for o = 1:nOutcomes
    training.(nameOutcomes{o}).outcome = outcomes.(nameOutcomes{o});
end

cd(pathLR), save('training','training')

cd(startpath)
end