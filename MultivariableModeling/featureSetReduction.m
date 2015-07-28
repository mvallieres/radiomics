function [fSet] = featureSetReduction(pathMINE,fSetNameType,outcome,setSize,nonTextStruct,textCells,textCellsName,paramAll,paramUsed,baseline,alpha,delta,nBoot,batchNum)
% -------------------------------------------------------------------------
% function [fSet] = buildFeatureSet(fSetName,outcome,setSize,nonTextStruct,textCells,paramAll,paramUsed,baseline,alpha,delta,nBoot)
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
% - pathMINE: Full path to the MINE.jar executable. The executable can be
%             downloaded at: <http://www.exploredata.net/>.
% - fSetNameType: String specifying the name of the type of feature set 
%                 (e.g., 'PET', 'SEPARATE', 'FUSED', etc.)
% - outcome: Column vector of size [nInst X 1] specifying the outcome status 
%            (1 or 0) for all instances.
% - setSize: Size of the output feature set (typically set to 25 in ref. [1]).
% - nonTextStruct: Structure data for non-texture features. This structure 
%                  is of the same format as the one saved as output to 
%                  computeAllNonTextureFeatures_STS.m, for example.                                
% - textCells: Cell vector of organized texture cells for all the
%              different scans. Format: textCells = {cell1, cell2, etc.}, 
%              where cell1, cell2, etc. are the files saved as output to
%              organizeSeparateTextures_STS.m or organizeFusedTextures_STS.m,
%              for example.
% - textCellsName: Cell vector corresponding to the name of the
%                  corresponding cells in textCells. 
%                  (e.g., textCellsName = {'PET','T1','T2FS'} for SEPARATE
%                  scans,  textCellsName = {'PET_T1','PET_T2FS'} for FUSED 
%                  scans, as defined in ref. [1])
% - paramAll: Cell vector incorporating all texture extraction parameters
%             tested in textCells. See EXAMPLE below for more details.
% - paramUsed: Vector of 1's and 0's to specify the degree of freedom on 
%              texture extraction parameters for the current experiment. 
%              For example, for an experiment where extraction parameters 
%              1, 2 and 4 in paramAll are allowed to vary, use
%              paramUsed = [1,1,0,1].
% - baseline: Vector of numerical values specifying the baseline texture 
%             extraction parameters for each entry in paramAll. See EXAMPLE
%             below for more details.
% - alpha: Numerical values specifying the coefficient of the first part of
%          the Gain equation, as defined in ref. [1].
% - delta: Numerical values specifying the coefficient of the second part 
%          of the Gain equation, as defined in ref. [1] (third part is set
%          to 0 in this function).
% - nBoot: Number of bootstrap samples to use.
% - batchNum: (optional input). If present, integer that specifies the
%             batch number for parallelization purposes
%
% See <https://github.com/mvallieres/radiomics/tree/master/STUDIES/STSstudy/Functions>
% to find computeAllNonTextureFeatures_STS.m, organizeSeparateTextures_STS.m
% and organizeFusedTextures_STS.m. See masterScript_STS.m for a complete 
% example of how to utilize the current function.
% -------------------------------------------------------------------------
% OUTPUTS:
% - fSET: Structure specifying the resulting feature set for the current
%         experiment.
%        --> fSET.Data: Array of size [nPatient X setSize], specifying the
%                       numerical data of the chosen features, where 'nInst'
%                       refers to the number of instances.
%        --> fSET.Info: Vector of cells with strings specifying information
%                       about the chosen features.
% -------------------------------------------------------------------------
% EXAMPLE:
% MRIinv_cell = {'NoInv','Inv'};
% MRIweight_mat = [1/4,1/3,1/2,2/3,3/4];
% R_mat = [1/2,2/3,1,3/2,2];
% scale_cell = {'pixelW',1,2,3,4,5};
% algo_cell = {'Equal','Lloyd'};
% Ng_mat = [8,16,32,64];
%
% FOR FUSED SCANS
% paramAll = {MRIinv_cell,MRIweight_mat,R_mat,scale_cell,algo_cell,Ng_mat};
% paramUsed = [1 1 0 0 0 1]; (example for a given experiment)
% baseline = [1 3 3 1 2 3];
%
% FOR SEPARATE SCANS
% paramAll = {R_mat,scale_cell,algo_cell,Ng_mat};
% paramUsed = [1 1 0 1]; (example for a given experiment)
% baseline = [3 1 2 3];
%
% NOTE: paramAll must always be of the same format and size for SEPARATE 
% and FUSED SCANS, with the same ordering of different extraction
% parameters.
% -------------------------------------------------------------------------
% AUTHOR(S): Martin Vallieres <mart.vallieres@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: May 2015
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


% RANDOM NUMBER GENERATOR SEED
if ~RandStream.getGlobalStream.Seed
    rng('shuffle')
    if nargin == 14
        % To avoid similar seeds when different batch are started with minimal time delay
        RandStream.setGlobalStream(RandStream('mt19937ar','seed',RandStream.getGlobalStream.Seed/(batchNum)^3))
    end
end


% INITIALIZATION
nInst = numel(outcome); % Number of patients
nScan = numel(textCellsName); % Number of scans included in that feature set experiment
fSet = struct; fSet.Data=zeros(nInst,setSize); fSet.Info=cell(setSize,1); % Initializing final feature set structure
if nargin == 14
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



% GETTING BOOTSTRAP SAMPLES TO BE USED FOR ALL EXPERIMENTS (using imbalance-adjusted resampling)
[bootSam,~] = buildBootSet(outcome,nBoot,'adjust');



% GETTING BOOTSTRAP RESULTS FOR THE FIRST PART OF THE GAIN EQUATION (ref. [1])
matCorr = zeros(nFeatures,nParam); % Matrix of Spearman correlation of each feature with the outcome (part 1 of Gain equation)
if sum(paramUsed(:)) % Need a transpose 
    for n = 1:nBoot
        dataBoot = matData(:,:,bootSam(:,n));
        outcomeBoot = outcome(bootSam(:,n));
        for i = 1:size(dataBoot,1)
            matCorr(i,:) = matCorr(i,:) + corr(squeeze(dataBoot(i,:,:))',outcomeBoot,'type','Spearman')';
        end
    end
else % No transpose
    for n = 1:nBoot
        dataBoot = matData(:,:,bootSam(:,n));
        outcomeBoot = outcome(bootSam(:,n));
        for i = 1:size(dataBoot,1)
            matCorr(i,:) = matCorr(i,:) + corr(squeeze(dataBoot(i,:,:)),outcomeBoot,'type','Spearman')';
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
    if numel(paramSize) == 4 % SEPARATE scans
        [ind1,ind2,ind3,ind4] = ind2sub(paramSize,col);
        indBase = [ind1,ind2,ind3,ind4];
        indFinal = zeros(1,4);
        for i = 1:4
            if isempty(param{i})
                indFinal(i) = baseline(i);
            else
                indFinal(i) = indBase(i);
            end
        end
        string = [cellNames{row,1},'(R=',num2str(paramAll{1}(indFinal(1)),'%.2f'),',Scale=',num2str(paramAll{2}{indFinal(2)}),',Quant.algo=',paramAll{3}{indFinal(3)},',Ng=',num2str(paramAll{4}(indFinal(4))),')--',cellNames{row,2}];
    elseif numel(paramSize) == 6 % FUSED scans
        [ind1,ind2,ind3,ind4,ind5,ind6] = ind2sub(paramSize,col);
        indBase = [ind1,ind2,ind3,ind4,ind5,ind6];
        indFinal = zeros(1,6);
        for i = 1:6
            if isempty(param{i})
                indFinal(i) = baseline(i);
            else
                indFinal(i) = indBase(i);
            end
        end
        string = [cellNames{row,1},'(',name,'Inv=',paramAll{1}{indFinal(1)},',',name,'w=',num2str(paramAll{2}(indFinal(2)),'%.2f'),',R=',num2str(paramAll{3}(indFinal(3)),'%.2f'),',Scale=',num2str(paramAll{4}{indFinal(4)}),',Quant.algo=',paramAll{5}{indFinal(5)},',Ng=',num2str(paramAll{6}(indFinal(6))),')--',cellNames{row,2}];
    end
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
        PICtest(:,f) = PICtest(:,f) + (1 - applyMIC(pathMINE,varBoot,dataBootAv,micName));
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
        if numel(paramSize) == 4 % SEPARATE scans
            [ind1,ind2,ind3,ind4] = ind2sub(paramSize,col);
            indBase = [ind1,ind2,ind3,ind4];
            indFinal = zeros(1,4);
            for i = 1:4
                if isempty(param{i})
                    indFinal(i) = baseline(i);
                else
                    indFinal(i) = indBase(i);
                end
            end
            string = [cellNames{row,1},'(R=',num2str(paramAll{1}(indFinal(1)),'%.2f'),',Scale=',num2str(paramAll{2}{indFinal(2)}),',Quant.algo=',paramAll{3}{indFinal(3)},',Ng=',num2str(paramAll{4}(indFinal(4))),')--',cellNames{row,2}];
        elseif numel(paramSize) == 6 % FUSED scans
            [ind1,ind2,ind3,ind4,ind5,ind6] = ind2sub(paramSize,col);
            indBase = [ind1,ind2,ind3,ind4,ind5,ind6];
            indFinal = zeros(1,6);
            for i = 1:6
                if isempty(param{i})
                    indFinal(i) = baseline(i);
                else
                    indFinal(i) = indBase(i);
                end
            end
            string = [cellNames{row,1},'(',name,'Inv=',paramAll{1}{indFinal(1)},',',name,'w=',num2str(paramAll{2}(indFinal(2)),'%.2f'),',R=',num2str(paramAll{3}(indFinal(3)),'%.2f'),',Scale=',num2str(paramAll{4}{indFinal(4)}),',Quant.algo=',paramAll{5}{indFinal(5)},',Ng=',num2str(paramAll{6}(indFinal(6))),')--',cellNames{row,2}];
        end
    else % This is a nonTexture feature
        string = cellNames{row,1};
    end
    fSet.Info{f+1} = string;
end

end