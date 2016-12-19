function [results] = predictionPerformanceEstimation_STS(X,Y,nBoot,imbalance,batchNum)
% -------------------------------------------------------------------------
% function [results] = predictionPerformanceEstimation(X,Y,nBoot,imbalance,batchNum)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function computes prediction performance estimation in terms of AUC, 
% sensitivity, specificity and accuracy according to the 0.632+ bootstrap 
% methodology. See ref. [1] for more details. This function uses logistic 
% regression utilities from DREES <http://www.cerr.info/drees>, and a fast 
% implementation of AUC calculation by Enric Junqué de Fortuny that is 
% available at: <http://www.mathworks.com/matlabcentral/fileexchange/41258-faster-roc-auc>
% -------------------------------------------------------------------------
% REFERENCE:
% [1] Vallieres, M. et al. (2015). A radiomics model from joint FDG-PET and 
%     MRI texture features for the prediction of lung metastases in soft-tissue 
%     sarcomas of the extremities. Physics in Medicine and Biology, 60(14), 
%     5471-5496. doi:10.1088/0031-9155/60/14/5471 
% [2] Vallieres, M. et al. (2015).
% -------------------------------------------------------------------------
% INPUTS:
% - X: Matrix of size [nInst X nFeat], specifying the numerical data of the 
%      features of the input features, where 'nInst' refers to the number 
%      of instances in X, and 'nFeat' to the number of features in X. 
%      Each column is a different feature.
% - Y: Column vector of size [nInst X 1] specifying the outcome status 
%      (1 or 0) for all instances.
% - nBoot: Number of bootstrap samples to use.
% - imbalance: String specifying the type of imbalance-adjustement strategy
%              employed. Either 'IABR' for imbalance-adjusted bootstrap
%              resampling (see ref.[1]), or 'IALR' for imbalance-adjusted
%              logistic regression (see ref.[2]).
% - batchNum: (optional input). If present, integer that specifies the
%             batch number for parallelization purposes
% -------------------------------------------------------------------------
% OUTPUTS:
% - results: Structure specifying the prediction performance estimation
%            results of the input multivariable models in terms of AUC, 
%            sensitivity, specificity and accuracy according to the 0.632+ 
%            bootstrap methodology and the ordinary bootstrap methods. 
%            Standard errors (SE) on a 95% condifence interval over all 
%            bootstrap samples are also provided.
% -------------------------------------------------------------------------
% AUTHOR(S): 
% - Martin Vallieres <mart.vallieres@gmail.com>
% - DREES development team <http://www.cerr.info/drees> (logistic regression)
% - Enric Junqué de Fortuny (fastAUC.cpp)
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: May 2015
% - Revision I: July 2015 (including imbalance-adjusted logistic regression) 
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
if ~RandStream.getGlobalStream.Seed
    rng('shuffle')
    if nargin == 5
        % To avoid similar seeds when different batch are started with minimal time delay
        RandStream.setGlobalStream(RandStream('mt19937ar','seed',RandStream.getGlobalStream.Seed/(batchNum)^3))
    end
end

% INITIALIZATION
top = 1-1/exp(1);
low = 1/exp(1);
bound = 1.96; % 95% confidence interval
results = struct;
AUC=zeros(1,nBoot);AUC632=zeros(1,nBoot);
SENS=zeros(1,nBoot);SENS632=zeros(1,nBoot);
SPEC=zeros(1,nBoot);SPEC632=zeros(1,nBoot);
ACCU=zeros(1,nBoot);ACCU632=zeros(1,nBoot);


% GETTING BOOTSTRAP SAMPLES FOR ALL EXPERIMENTS + LOGISTIC REGRESSION TYPE
if strcmp(imbalance,'IABR')
    logisticRegression = @(x,y) applyStandardLR(x,y);
    [trainSets,testSets] = buildBootSet(Y,nBoot,'adjust'); % 'trainSets' is a matrix, 'testSets' is a cell 
elseif strcmp(imbalance,'IALR')
    logisticRegression = @(x,y) applyEnsembleLR(x,y);
    [trainSets,testSets] = buildBootSet(Y,nBoot); % 'trainSets' is a matrix, 'testSets' is a cell 
end



% COMPUTING RESULTS FOR THE WHOLE DATASET
[Xtrain] = normalizeZeroOne(X);
[coeff] = logisticRegression(Xtrain,Y);
[resp] = responseLR(Xtrain,coeff);
[aucData,sensData,specData,accuData] = calcPerformMetrics(resp,Y,0);


% COMPUTING RESULTS FOR BOOTSTRAP SAMPLES
for n = 1:nBoot
    Xtrain = X(trainSets(:,n),:); Xtest = X(testSets{n},:); Ytrain = Y(trainSets(:,n),1); Ytest = Y(testSets{n},1);
    [Xtrain,Xtest] = normalizeZeroOne(Xtrain,Xtest);
    [coeff] = logisticRegression(Xtrain,Ytrain);
    [resp] = responseLR(Xtest,coeff);
    [aucBoot,sensBoot,specBoot,accuBoot] = calcPerformMetrics(resp,Ytest,0);
    AUC(n) = aucBoot; SENS(n) = sensBoot; SPEC(n) = specBoot; ACCU(n) = accuBoot;
    
    % For AUC
    alpha = top/(1-low*(aucData-aucBoot)/(aucData-0.5+eps));
    if alpha > 1
       alpha = 1;
    elseif alpha < top
       alpha = top;
    end
    if aucBoot < 0.5
       aucBoot = 0.5;
    end
    AUC632(n) = (1-alpha)*aucData+alpha*aucBoot;
    
    % For Sens
    alpha = top/(1-low*(sensData-sensBoot)/(sensData+eps));
    if alpha < top
        alpha = top;
    end
    SENS632(n) = (1-alpha)*sensData+alpha*sensBoot;
    
    % For Spec
    alpha = top/(1-low*(specData-specBoot)/(specData+eps));
    if alpha < top
        alpha = top;
    end
    SPEC632(n) = (1-alpha)*specData+alpha*specBoot;
    
    % For Accu
    alpha = top/(1-low*(accuData-accuBoot)/(accuData+eps));
    if alpha < top
        alpha = top;
    end
    ACCU632(n) = (1-alpha)*accuData+alpha*accuBoot;
end
results.AUC = mean(AUC); results.AUC632 = mean(AUC632);
results.Sensitivity = mean(SENS); results.Sensitivity632 = mean(SENS632);
results.Specificity = mean(SPEC); results.Specificity632 = mean(SPEC632);
results.Accuracy = mean(ACCU); results.Accuracy632 = mean(ACCU632);
results.SE_AUC = bound*std(AUC)/sqrt(nBoot); results.SE_AUC632 = bound*std(AUC632)/sqrt(nBoot);
results.SE_Sensitivity = bound*std(SENS)/sqrt(nBoot); results.SE_Sensitivity632 = bound*std(SENS632)/sqrt(nBoot);
results.SE_Specificity = bound*std(SPEC)/sqrt(nBoot); results.SE_Specificity632 = bound*std(SPEC632)/sqrt(nBoot);
results.SE_Accuracy = bound*std(ACCU)/sqrt(nBoot); results.SE_Accuracy632 = bound*std(ACCU632)/sqrt(nBoot);
    
end