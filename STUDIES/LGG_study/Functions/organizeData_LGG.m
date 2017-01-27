function organizeData_LGG(pathTCGA,pathSTUDY,patientID,scans,outcomesTCGA,timeToEventTCGA,textType,textName,paramAll)
% -------------------------------------------------------------------------
% function organizeData_LGG(pathTCGA,pathSTUDY,patientID,scans,outcomesTCGA,timeToEventTCGA,textType,textName,paramAll);
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function organizes radiomic data for available patients for each
% outcome studied in the LGG study.
% -------------------------------------------------------------------------
% INPUTS:
% 1. pathTCGA: Full path to the 'TCGA_DATA' of the LGG study. Defined in
%              masterScript_LGG.m. 
%              --> Ex: '/myProject/WORKSPACE/TCGA_DATA'
% 2. pathSTUDY: Full path to where the organized data will be saved.
%               Defined in masterScript_LGG.m
%               --> Ex: '/myProject/WORKSPACE/STUDY_DATA'
% 3. patientID: Cell of strings of size {nPatient X 1} specifying the IDs
%               of the LGG patients used in this study.
%               --> Ex: {108 X 1} cell of strings
% 4. scans: Cell of strings specifying the type of MRI scans used in this
%           study.
%           --> Ex: {'T1W','T1CE','T2W','T2F'}
% 5. outcomesTCGA: Structure specifying the status (1 or 0) for different
%                  outcomes in LGG cancer. Contains: outcomes.nonIDH1, 
%                  outcomes.IDHcodel, outcomes.progression and
%                  outcomes.lowGrade.
%                  --> Ex: Structure containing 4 [nPatient X 1] vectors
% 6. timeToEventTCGA: Structure specifying the time to events for different
%                     outcomes in LGG cancer. 
%                     Contains: timeToEvent.progressionFreeTime 
%                     --> Ex: Structure containing 1 [nPatient X 1] vector
% 7. textType: Cell of strings specifying the texture types used in this study.
%              --> Defined in masterScript_LGG.m
% 8. textName: Cell of strings specifying the texture names used in this study.
%              --> Defined in masterScript_LGG.m
% 9. paramAll: Cell specifying all the different texture extraction
%              parameters used in this study.
%              --> Defined in masterScript_LGG.m
% -------------------------------------------------------------------------
% OUTPUTS: Organized 'outcomes', 'timeToEvent', 'lesionNum' and texture and 
%          non-texture data in 'pathSTUDY'. 
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
cd(pathTCGA)


% INITIALIZATION
load('nonTextures') % Variable 'nonText' now in MATLAB workspace
cd('TEXTURES'), pathText = pwd;
nScans = numel(scans);
nameOutcomes = fieldnames(outcomesTCGA); nOutcomes = numel(nameOutcomes);
nonTextName = fieldnames(nonText); nNonText = numel(nonTextName);
nPatient = numel(patientID);
nParamType = max(size(paramAll));
nParam = zeros(1,nParamType);
for i = 1:nParamType
    nParam(i) = numel(paramAll{i});
end
for o = 1:nOutcomes
    lesionNum.(nameOutcomes{o}) = [];
end
excludeScans = 0;
excludeOutcomes = struct;
for o = 1:nOutcomes
    excludeOutcomes.(nameOutcomes{o}) = 0;
end


% ANALYZING WHICH PATIENTS HAVE COMPLETE DATA
cd(pathTCGA), fid = fopen('exclusion_report.txt','w'); cd(pathText)
for p = 1:nPatient
    ok = 1;
    for s = 1:nScans
        try
            textures = load([patientID{p},'_',scans{s}]); textures = struct2cell(textures); textures = textures{1};
        catch
            ok = 0;
        end
    end
    for n = 1:nNonText
        if isnan(nonText.(nonTextName{n})(p))
            ok = 0;
        end
    end
    if ok
        for o = 1:nOutcomes
            if ~isnan(outcomesTCGA.(nameOutcomes{o})(p))
                lesionNum.(nameOutcomes{o}) = [lesionNum.(nameOutcomes{o});p];
            else
                excludeOutcomes.(nameOutcomes{o}) = excludeOutcomes.(nameOutcomes{o}) + 1;
            end
        end
    else
        excludeScans = excludeScans + 1;
    end
end
fprintf(fid,'\n\n -------------------- TOTAL EXLUSION --------------------\n');
fprintf(fid,'--> INITIAL NUMBER OF PATIENTS: %u\n',nPatient);
fprintf(fid,'\n');
fprintf(fid,'--> OUT OF %u TOTAL PATIENTS, MISSING EITHER T1W, T2W, T1CE OR T2F: %u patients\n',nPatient,excludeScans);
fprintf(fid,'--> REMAINING PATIENTS: %u\n',nPatient-excludeScans);
fprintf(fid,'\n');
for o = 1:nOutcomes
    fprintf(fid,'--> OUT OF %u REMAINING PATIENTS, ALSO MISSING OUTCOME ''%s'': %u patients\n',nPatient-excludeScans,nameOutcomes{o},excludeOutcomes.(nameOutcomes{o}));
end
fprintf(fid,'\n');
for o = 1:nOutcomes
    fprintf(fid,'--> TOTAL REMAINING FOR OUTCOME ''%s'': %u patients\n',nameOutcomes{o},nPatient - excludeScans - excludeOutcomes.(nameOutcomes{o}));
end
fprintf(fid,'-------------------------------------------------------------\n\n');
fclose(fid);


% ARRANGING DATA

% Outcomes
outcomes = struct; timeToEvent = struct;
for o = 1:nOutcomes
    outcomes.(nameOutcomes{o}) = outcomesTCGA.(nameOutcomes{o})(lesionNum.(nameOutcomes{o}));
end
timeToEvent.progressionFreeTime = timeToEventTCGA.progressionFreeTime(lesionNum.progression);
cd(pathSTUDY), save('outcomes','outcomes'), save('timeToEvent','timeToEvent'), save('lesionNum','lesionNum')

% Non-Textures
tempNonText = nonText;
for o = 1:nOutcomes
    clear nonText
    for n = 1:nNonText
        nonText.(nonTextName{n}).Data = tempNonText.(nonTextName{n})(lesionNum.(nameOutcomes{o}));
        [nonText.(nonTextName{n}).Spearman.rs,nonText.(nonTextName{n}).Spearman.p] = corr(nonText.(nonTextName{n}).Data,outcomes.(nameOutcomes{o}),'type','Spearman','rows','pairwise');
    end
    cd(pathSTUDY), save(['nonText_',nameOutcomes{o}],'nonText')
end

% Textures
scale_mat = paramAll{1}; algo_cell = paramAll{2}; Ng_mat = paramAll{3};
for s = 1:nScans
    for o = 1:nOutcomes
        cd(pathText)
        nPatient = numel(lesionNum.(nameOutcomes{o}));
        text = cell(nParam);
        for n = 1:prod(nParam)
            for t = 1:numel(textType)
                for f = 1:numel(textName{t})
                    text{n}.(textType{t}).(textName{t}{f}).Data = zeros(nPatient,1);
                end
            end
        end
        for p = 1:nPatient
            textures = load([patientID{lesionNum.(nameOutcomes{o})(p)},'_',scans{s}]); textures = struct2cell(textures); textures = textures{1};
            exp = 0;
            for ss = 1:numel(scale_mat)
                for a = 1:numel(algo_cell)
                    for n = 1:numel(Ng_mat)
                        exp = exp + 1;
                        for t = 1:numel(textType)
                            for f = 1:numel(textName{t})
                                text{ss,a,n}.(textType{t}).(textName{t}{f}).Data(p) = textures.(['Experiment',num2str(exp)]).(textType{t}).(textName{t}{f});
                            end
                        end
                    end
                end
            end
        end
        for ss = 1:numel(scale_mat)
            for a = 1:numel(algo_cell)
                for n = 1:numel(Ng_mat)
                    for t = 1:numel(textType)
                        for f = 1:numel(textName{t})
                            [text{ss,a,n}.(textType{t}).(textName{t}{f}).Spearman.rs,text{ss,a,n}.(textType{t}).(textName{t}{f}).Spearman.p] = corr(text{s,a,n}.(textType{t}).(textName{t}{f}).Data,outcomes.(nameOutcomes{o}),'type','Spearman','rows','pairwise');
                        end
                    end
                end
            end
        end
        cd(pathSTUDY), save(['text_',nameOutcomes{o},'_',scans{s}],'text')
    end
end

% Vasari
cd(pathTCGA)
vasariInit = load('vasari'); vasariInit = struct2cell(vasariInit); vasariInit = vasariInit{1}; 
cd(pathSTUDY)
for o = 1:nOutcomes
    vasari = vasariInit(lesionNum.(nameOutcomes{o}),:);
    save(['vasari_',nameOutcomes{o}],'vasari')
end

cd(startpath)
end