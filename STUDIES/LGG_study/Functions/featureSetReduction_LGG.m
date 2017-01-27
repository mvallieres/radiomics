function [fSet] = featureSetReduction_LGG(pathMINE,outcome,setSize,nonTextStruct,textCells,textCellsName,paramAll,paramUsed,baseline,alpha,delta,nBoot,seed,batchNum)
% -------------------------------------------------------------------------
% function [fSet] = featureSetReduction_LGG(pathMINE,outcome,setSize,nonTextStruct,textCells,textCellsName,paramAll,paramUsed,baseline,alpha,delta,nBoot,seed,batchNum)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function computes feature set reduction for an experiment with a 
% specific degree of freedom on texture extraction parameters as defined by 
% 'paramUsed', according to the methodology described in ref. [1].
% -------------------------------------------------------------------------
% REFERENCE:
% [1] Vallieres, M. et al. (2015). A radiomics model from joint FDG-PET and 
%     MRI texture features for the prediction of lung metastases in soft-tissue 
%     sarcomas of the extremities. Physics in Medicine and Biology, 60(14), 
%     5471-5496. doi:10.1088/0031-9155/60/14/5471
% -------------------------------------------------------------------------
% INPUTS:
% 1. pathMINE: Full path to the MINE.jar executable. The executable can be
%              downloaded at: <http://www.exploredata.net/>.
% 2. outcome: Column vector of size [nInst X 1] specifying the outcome status 
%             (1 or 0) for all instances.
% 3. setSize: Size of the output feature set (typically set to 25 in ref. [1]).
% 4. nonTextStruct: Structure data for non-texture features. This structure 
%                   is of the same format as the one saved as output to 
%                   computeAllNonTextureFeatures_LGG.m, for example.                                
% 5. textCells: Cell vector of organized texture cells for all the
%               different scans. Format: textCells = {cell1, cell2, etc.}, 
%               where cell1, cell2, etc. are the files saved as output to
%               organizeData_LGG.m
% 6. textCellsName: Cell of strings corresponding to the name of the
%                   corresponding cells in textCells. 
% 7. paramAll: Cell vector incorporating all texture extraction parameters
%              tested in textCells. See EXAMPLE below for more details.
% 8. paramUsed: Vector of 1's and 0's to specify the degree of freedom on 
%               texture extraction parameters for the current experiment. 
%               For example, for an experiment where extraction parameters 
%               1, 2 and 3 in paramAll are allowed to vary, use
%               paramUsed = [1,1,1].
% 9. baseline: Vector of numerical values specifying the baseline texture 
%              extraction parameters for each entry in paramAll. See EXAMPLE
%              below for more details.
% 10. alpha: Numerical values specifying the coefficient of the first part of
%            the Gain equation, as defined in ref. [1].
% 11. delta: Numerical values specifying the coefficient of the second part 
%            of the Gain equation, as defined in ref. [1] (third part is set
%            to 0 in this function).
% 12. nBoot: Number of bootstrap samples to use.
% 13. batchNum: (optional input). If present, integer that specifies the
%               batch number for parallelization purposes
%
% See masterScript_LGG.m for a complete example of how to utilize the 
% current function.
% -------------------------------------------------------------------------
% OUTPUTS:
% 1. fSET: Structure specifying the resulting feature set for the current
%          experiment.
%          --> fSET.Data: Array of size [nPatient X setSize], specifying the
%                         numerical data of the chosen features, where 'nInst'
%                         refers to the number of instances.
%          --> fSET.Info: Vector of cells with strings specifying information
%                         about the chosen features.
% -------------------------------------------------------------------------
% EXAMPLE:
% scale_mat = [1,2,3,4,5];
% algo_cell = {'Equal','Uniform'};
% Ng_mat = [8,16,32,64];
%
% paramAll = {scale_mat,algo_cell,Ng_mat};
% paramUsed = [1 1 1]; (example for a given experiment)
% baseline = [1,2,3];
% -------------------------------------------------------------------------
% AUTHOR(S): Martin Vallieres <mart.vallieres@gmail.com>
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


% RANDOM NUMBER GENERATOR SEED
rng(seed);


% INITIALIZATION
nInst = numel(outcome); % Number of patients
nScan = numel(textCellsName); % Number of scans included in that feature set experiment
fSet = struct; fSet.Data=zeros(nInst,setSize); fSet.Info=cell(setSize,1); % Initializing final feature set structure
if nargin == 15
    micName = ['batch',num2str(batchNum)];
else
    micName = 'master';
end
if ~isempty(strfind(textCellsName{1},'CT'))
    name = 'CT';
else
    name = 'MR';
end

% Initializing names
nonTextName = fieldnames(nonTextStruct); nNonText = numel(nonTextName);
textTypeName = fieldnames(textCells{1}{1}); nTextType = numel(textTypeName);
nameF = cell(1,nTextType);
nText = 0;
for t = 1: nTextType
    nameF{t} = fieldnames(textCells{1}{1}.(textTypeName{t}));
    nText = nText + numel(nameF{t});
end
textName = cell(nText,1);
count = 1;
for t = 1:nTextType
    for f = 1:numel(nameF{t})
        textName{count} = [textTypeName{t},'-',nameF{t}{f}];
        count = count +1;
    end
end
nFeatures = nNonText + nScan*nText;
cellNames = cell(nFeatures,2);
count = 1;
for i = 1:nScan
    for j = 1:nText
        cellNames{count,1} = textCellsName{i};
        cellNames{count,2} = textName{j};
        count = count + 1;
    end
end
for i = 1:nNonText
    cellNames{count,1} = nonTextName{i};
    count = count + 1;
end

% Extraction parameters used in that experiment
nParamType = numel(paramUsed);
param = cell(1,nParamType); 
paramSize = ones(1,nParamType);
nParam = 1;
for i = 1:nParamType
    if paramUsed(i)
        temp = length(paramAll{i});
        nParam = nParam * temp;
        paramSize(i) = temp;
        param{i} = paramAll{i};
    end
end



% FILLING UP WORKING DATA MATRIX
matData = zeros(nFeatures,nParam,nInst); % Matrix containing all the data vectors for each feature (used for part 2 of Gain equation)
count = 1;
for i = 1:nScan
    textCell = textCells{i};
    textCell = removeParam(textCell,paramUsed,baseline);
    for t = 1:nTextType
        for f = 1:numel(nameF{t})
            for p = 1:numel(textCell)
                matData(count,p,:) = textCell{p}.(textTypeName{t}).(nameF{t}{f}).Data;
            end
            count = count + 1;
        end
    end
end
for i = 1:nNonText
    matData(nText*nScan + i,1,:) = nonTextStruct.(nonTextName{i}).Data;
    for j = 2:nParam
        matData(nText*nScan + i,j,:) = matData(nText*nScan + i,1,:);
    end
end


% IMPORTANT: IF NAN, IT WOULD MEAN SOMETHING WENT WRONG IN SOME TEXTURE CALCULATIONS
matData(isnan(matData(:))) = 0;


% GETTING BOOTSTRAP SAMPLES TO BE USED FOR ALL EXPERIMENTS (using imbalance-adjusted resampling)
[bootSam,~] = buildBootSet(outcome,nBoot,'adjust');



% GETTING BOOTSTRAP RESULTS FOR THE FIRST PART OF THE GAIN EQUATION (ref. [1])
matCorr = zeros(nFeatures,nParam); % Matrix of Spearman correlation of each feature with the outcome (part 1 of Gain equation)
if sum(paramUsed(:)) % Need a transpose 
    for n = 1:nBoot
        dataBoot = matData(:,:,bootSam(:,n));
        outcomeBoot = outcome(bootSam(:,n));
        for i = 1:size(dataBoot,1)
            matCorr(i,:) = matCorr(i,:) + corr(squeeze(dataBoot(i,:,:))',outcomeBoot,'type','Spearman','rows','pairwise')';
        end
    end
else % No transpose
    for n = 1:nBoot
        dataBoot = matData(:,:,bootSam(:,n));
        outcomeBoot = outcome(bootSam(:,n));
        for i = 1:size(dataBoot,1)
            matCorr(i,:) = matCorr(i,:) + corr(squeeze(dataBoot(i,:,:)),outcomeBoot,'type','Spearman','rows','pairwise')';
        end
    end
end
matCorr = abs(matCorr./nBoot);



% CHOOSING FIRST FEATURE(depends only on Spearman's correlation)
indChosen = zeros(setSize-1,1);

% Choosing feature
[~,index] = max(matCorr(:));
[row,col] = ind2sub([nFeatures,nParam],index);
indChosen(1) = row;
fSet.Data(:,1) = squeeze(matData(row,col,:));

% Obtaining feature name
if row <= nScan*nText
    [ind1,ind2,ind3] = ind2sub(paramSize,col);
    indBase = [ind1,ind2,ind3];
    indFinal = zeros(1,3);
    for i = 1:3
        if isempty(param{i})
            indFinal(i) = baseline(i);
        else
            indFinal(i) = indBase(i);
        end
    end
    string = [cellNames{row,1},'(R=',num2str(1,'%.2f'),',Scale=',num2str(paramAll{1}(indFinal(1))),',Quant.algo=',paramAll{2}{indFinal(2)},',Ng=',num2str(paramAll{3}(indFinal(3))),')--',cellNames{row,2}];
else % This is a nonTexture feature
    string = cellNames{row,1};
end
fSet.Info{1} = string;



% COMPUTING FOR OTHER FEATURES
PICtest = zeros(nFeatures,setSize-1);
for f = 1:setSize-1
    
    % Computing PIC for all bootstrap samples
    for n = 1:nBoot
        varBoot = fSet.Data(bootSam(:,n),f);
        dataBoot = matData(:,:,bootSam(:,n));
        dataBootAv = squeeze(mean(dataBoot,2))';
        try
            PICtest(:,f) = PICtest(:,f) + (1 - applyMIC(pathMINE,varBoot,dataBootAv,micName));
        catch % SOLVE THAT PROBLEM
            PICtest(:,f) = PICtest(:,f) + zeros(nFeatures,1);
        end
    end
    PICtest(:,f) = PICtest(:,f)./nBoot;
    PICtemp = zeros(nFeatures,1);
    for k = 1:f
        PICtemp = PICtemp + 2*(f-k+1)/(f*(f+1)).*PICtest(:,k);
    end
    PIC = repmat(PICtemp,1,nParam);
    
    % Choosing feature
    Gain = alpha.*matCorr + delta.*PIC;
    [~,index] = sort(Gain(:),'descend'); best = 1;
    [row,col] = ind2sub([nFeatures,nParam],index(best));
    while ~isempty(find(indChosen==row))
        best = best + 1;
        [row,col] = ind2sub([nFeatures,nParam],index(best));
    end
    indChosen(f+1) = row;
    fSet.Data(:,f+1) = squeeze(matData(row,col,:));
    
    % Obtaining feature name
    if row <= nScan*nText
        [ind1,ind2,ind3] = ind2sub(paramSize,col);
        indBase = [ind1,ind2,ind3];
        indFinal = zeros(1,3);
        for i = 1:3
            if isempty(param{i})
                indFinal(i) = baseline(i);
            else
                indFinal(i) = indBase(i);
            end
        end
        string = [cellNames{row,1},'(R=',num2str(1,'%.2f'),',Scale=',num2str(paramAll{1}(indFinal(1))),',Quant.algo=',paramAll{2}{indFinal(2)},',Ng=',num2str(paramAll{3}(indFinal(3))),')--',cellNames{row,2}];
    else % This is a nonTexture feature
        string = cellNames{row,1};
    end
    fSet.Info{f+1} = string;
end

end