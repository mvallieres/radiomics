function [models] = featureSelectionRF_HN(X,categories,Y,maxOrder,nBoot,Info,imbalance,seed)
% -------------------------------------------------------------------------
% function [models] = featureSelection(X,Y,maxOrder,nBoot,Info,imbalance,batchNum)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function computes feature set selection according to the 0.632+ 
% bootstrap methodology for an input matrix of features and and input 
% outcome vector, and for multiple model orders as defined by the user. 
% See ref. [1] for more details. This function uses logistic regression 
% utilities from DREES <http://www.cerr.info/drees>, and a fast 
% implementation of AUC calculation by Enric Junqué de Fortuny that is
% available at: <http://www.mathworks.com/matlabcentral/fileexchange/41258-faster-roc-auc>
% 
% NOTE: This function now maximizes 0.5*AUC + 0.5*(1-abs(SENSITIVITY-SPECIFICITY))
% -------------------------------------------------------------------------
% REFERENCE:
% [1] Vallieres, M. et al. (2015). A radiomics model from joint FDG-PET and 
%     MRI texture features for the prediction of lung metastases in soft-tissue 
%     sarcomas of the extremities. Physics in Medicine and Biology, 60(14), 
%     5471-5496. doi:10.1088/0031-9155/60/14/5471
% -------------------------------------------------------------------------
% INPUTS:
% - X: Matrix of size [nInst X nFeat], specifying the numerical data of the 
%      features of the input features, where 'nInst' refers to the number 
%      of instances in X, and 'nFeat' to the number of features in X. 
%      Each column is a different feature.
% - Y: Column vector of size [nInst X 1] specifying the outcome status 
%      (1 or 0) for all instances.
% - maxOrder: Integer specifying the maximal model order to construct.
% - nBoot: Number of bootstrap samples to use.
% - Info: Cell of size [nFeat X 1] of strings specifying the name of each 
%         feature in 'X'.
% - imbalance: String specifying the type of imbalance-adjustement strategy
%              employed. Either 'IABR' for imbalance-adjusted bootstrap
%              resampling (see ref.[1]), or 'IALR' for imbalance-adjusted
%              logistic regression.
% - batchNum: (optional input). If present, integer that specifies the
%             batch number for parallelization purposes.
% -------------------------------------------------------------------------
% OUTPUTS:
% - models: Structure specifying the resulting multivariable models for the
%           input feature set in 'data', for each model order. Example for
%           order 4:
%        --> models.Order4.Data: Matrix of size [nInst X order], specifying
%                                the selected features, in order of
%                                selection.
%        --> models.Order4.Name: Cell specifying the names of the selected
%                                features, in order of selection.
% -------------------------------------------------------------------------
% AUTHOR(S): 
% - Martin Vallieres <mart.vallieres@gmail.com>
% - DREES development team <http://www.cerr.info/drees> (logistic regression)
% - Enric Junqué de Fortuny (fastAUC.cpp)
% -------------------------------------------------------------------------
% HISTORY:
% - Creation - May 2015
% - Revision I - July 2015: including imbalance-adjusted logistic regression 
% - Revision II - July 2015: maximizing 0.5*AUC + 0.5*(1-abs(SENSITIVITY-SPECIFICITY))
%--------------------------------------------------------------------------
% STATEMENT:
% This file is part of <https://github.com/mvallieres/radiomics/>, 
% a package providing MATLAB programming tools for radiomics analysis.
% --> Copyright (C) 2015  Martin Vallieres
% --> Copyright 2010, Joseph O. Deasy, on behalf of the DREES development team.
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
%
%    _______________________________________________________________
%
% --> Copyright (c) 2013, Enric Junqué de Fortuny
%     All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%      
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.
% -------------------------------------------------------------------------


% RANDOM NUMBER GENERATOR SEED
rng(seed);


% INITIALIZATION
top = 1-1/exp(1);
low = 1/exp(1);
nFeat = numel(X);
models = struct; % Final model structure
if nFeat < maxOrder
    maxOrder = nFeat;
end
nBoot = 100; nSplit = 10; testSplit = 1/3;


rng(seed), [trainSets,testSets] = getSplits(Y,nSplit,testSplit);


% FORWARD FEATURE SELECTION (for all different starters)
modelMat = zeros(nFeat,maxOrder);
metricMat = zeros(nFeat,maxOrder);
for i = 1:nFeat 
     indLeft = 1:nFeat;
     % Order 1
     indLeft(i) = [];
     modelMat(i,1) = i;   
  
     aucTemp = 0; sensTemp = 0; specTemp = 0;
     for n = 1:nSplit
         
         % Random forest
         nCol = 1; Ttrain = table; Ttest = table;
         for f = 1:nCol
             Ttrain = [Ttrain,table(X{i}(trainSets(:,n)),'VariableNames',{['var',num2str(f)]})];
             Ttest = [Ttest,table(X{i}(testSets(:,n)),'VariableNames',{['var',num2str(f)]})];
         end
         cat = categories(i);
         Ytrain = Y(trainSets(:,n),1); Ytest = Y(testSets(:,n),1);
         rng(seed), [RF] = trainRF_table(Ttrain,Ytrain,cat,nBoot);
         [prob] = predictRF(Ttest,RF);
         [aucSplit,sensSplit,specSplit,~] = calcPerformMetrics(prob,Ytest,0.5);

         aucTemp = aucTemp + aucSplit;
         sensTemp = sensTemp + sensSplit;
         specTemp = specTemp + specSplit;  
     end
     aucTemp = aucTemp/nSplit; sensTemp = sensTemp/nSplit; specTemp = specTemp/nSplit;
     metricMat(i,1) = 0.5*aucTemp + 0.5*(1-abs(sensTemp-specTemp));
     
     % Going for orders 2 to maxOrder
     for j = 2:maxOrder
         maxMetric = 0;
         for k = 1:(nFeat-j+1)
             indexModel = [modelMat(i,1:(j-1)),indLeft(k)];

             aucTemp = 0; sensTemp = 0; specTemp = 0;
             for n=1:nSplit
                 
                 % Random forest
                 nCol = numel(indexModel); Ttrain = table; Ttest = table;
                 for f = 1:nCol
                     Ttrain = [Ttrain,table(X{indexModel(f)}(trainSets(:,n)),'VariableNames',{['var',num2str(f)]})];
                     Ttest = [Ttest,table(X{indexModel(f)}(testSets(:,n)),'VariableNames',{['var',num2str(f)]})];
                 end
                 cat = categories(indexModel);
                 Ytrain = Y(trainSets(:,n),1); Ytest = Y(testSets(:,n),1);
                 rng(seed), [RF] = trainRF_table(Ttrain,Ytrain,cat,nBoot);
                 [prob] = predictRF(Ttest,RF);
                 [aucSplit,sensSplit,specSplit,~] = calcPerformMetrics(prob,Ytest,0.5); 
         
                 aucTemp = aucTemp + aucSplit;
                 sensTemp = sensTemp + sensSplit;
                 specTemp = specTemp + specSplit; 
             end
             aucTemp = aucTemp/nSplit; sensTemp = sensTemp/nSplit; specTemp = specTemp/nSplit;
             metricTemp = 0.5*aucTemp + 0.5*(1-abs(sensTemp-specTemp));
             if metricTemp >= maxMetric
                 maxMetric = metricTemp;
                 index = indLeft(k);
             end
         end
         modelMat(i,j) = index;
         metricMat(i,j) = maxMetric;
         indLeft(find(indLeft==index)) = [];
     end
end


% OBTAINING MAXIMUM RESULTS FOR EVERY MODEL ORDER (maximum from all different starters)
[~,indMax] = max(metricMat);
for i = 1:maxOrder
    nameOrder = ['Order',num2str(i)];
    models.(nameOrder).Data = X(modelMat(indMax(i),1:i));
    models.(nameOrder).Categories = categories(modelMat(indMax(i),1:i));
    models.(nameOrder).Name = cell(i,1);
    for j = 1:i
        models.(nameOrder).Name{j} = Info{modelMat(indMax(i),j)};
    end
end

end