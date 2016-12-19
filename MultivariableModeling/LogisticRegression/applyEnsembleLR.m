function [coeff] = applyEnsembleLR(X,Y,seed)
% -------------------------------------------------------------------------
% function [coeff] = applyEnsembleLR(X,Y,seed)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function computes imbalance-adjusted logistic regression (IALR) by
% creating an ensemble of balanced classifiers as shown in ref.[1].
% Basically, the data is partioned into a collection of balanced subsets,
% where each subset consists of all the positive instances and an equal
% number of negative instances drawn from a random permutation . For 
% example, let a data set containing 12 positive instances, and 145 
% negative instances (total: 157). The ensemble classifier would consists 
% of 12 subsets each containing all positive instances, where 12 subsets 
% would contain 12 negative instances, and 1 subset would contain 13 
% negative instances. A logistic regression classifier is then built for 
% each subset of the data, and the final coefficients are calculated as the 
% mean of the corresponding coefficients of each classifier.  This function 
% uses logistic regression utilities from DREES <http://www.cerr.info/drees>.
% -------------------------------------------------------------------------
% REFERENCE:
% [1] Schiller, T.W. et al. (2015). Modeling radiation-induced lung injury
%     risk with an ensemble of support vector machines. Neurocomputing, 73 
%     (10-12), 1861-1867.
% -------------------------------------------------------------------------
% INPUTS:
% - X: Matrix of size [nInst,nFeat], specifying the numerical data of the 
%      features of the input features, where 'nInst' refers to the number 
%      of instances in X, and 'nFeat' to the number of features in X. 
%      Each column is a different feature.
% - Y: Column vector of size [nInst,1] specifying the outcome status 
%      (1 or 0) for all instances.
% - seed: (optional). Random generator seed for reproducibility of
%         experiments.
% -------------------------------------------------------------------------
% OUTPUTS:
% - coeff: Column vector of size [nCoeff+1,1] specifying the final 
%          logistic regression coefficients. Last entry specify the offset
%          of the model.
% -------------------------------------------------------------------------
% AUTHOR(S): 
% - Martin Vallieres <mart.vallieres@gmail.com>
% - DREES development team <http://www.cerr.info/drees> (logistic regression)
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: December 2016
%--------------------------------------------------------------------------
% STATEMENT:
% This file is part of <https://github.com/mvallieres/radiomics/>, 
% a package providing MATLAB programming tools for radiomics analysis.
% --> Copyright (C) 2015-2016  Martin Vallieres
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
% -------------------------------------------------------------------------


% INITIALIZATION
nNeg = sum(~Y); % Number of negative instances
nPos = sum(Y); % Number of positive instances
indPos = find(Y); % Indexes of positive instances
indNeg = find(~Y); % Indexes of negative instances
nSub = round(nNeg/nPos); % Number of subsets


% RANDOM NUMBER GENERATOR SEED
if nargin == 3
    rng(seed)
else
    if ~RandStream.getGlobalStream.Seed
        rng('shuffle')
    end
end

% FIND THE NUMBER OF NEGATIVE INSTANCES IN EACH SUBSET
nNegS = nNeg / nSub; nSup = ceil(nNegS); nInf = floor(nNegS);
if nSup ~= nInf
    nSubInf = nSub - 1; nSubSup = 1; total = nSubInf*nInf + nSubSup*nSup;
    while total ~= nNeg
        nSubInf = nSubInf - 1; nSubSup = nSubSup + 1;
        total = nSubInf*nInf + nSubSup*nSup;
    end
    nNegS = [repmat(nInf,[1,nSubInf]),repmat(nSup,[1,nSubSup])];
else % The number of negative instances in all partitions will be the same
    nNegS = repmat(nSup,[1,nSub]);
end


% FINDING THE INDEXES OF NEGATIVE INSTANCES IN EACH SUBSET
indNegSub = cell(1,nSub);
for i = 1:nSub-1
    indNegSub{i} = zeros(nNegS(i),1);
    indTemp = ceil(numel(indNeg)*rand(nNegS(i),1));
    indTemp = unique(indTemp);
    total = numel(indTemp);
    while total ~= nNegS(i)
        indMore = ceil(numel(indNeg)*rand(nNegS(i)-total,1));
        indTemp = [indTemp;unique(indMore)];
        indTemp = unique(indTemp);
        total = numel(indTemp);
    end
    indNegSub{i}(:) = indNeg(indTemp);
    indNeg(indTemp) = [];
end
indNegSub{end} = indNeg;


% COMPUTING LOGISTIC REGRESSION FOR EACH SUBSET AND AVERAGING
coeff = zeros(size(X,2)+1,1);
for i = 1:nSub
    ind = [indNegSub{i};indPos];
    Xtemp = X(ind,:);
    Ytemp = Y(ind);
    [Xtemp,Ytemp] = shufflePartition(Xtemp,Ytemp);
    try
        [~,temp,~] = drxlr_apply_logistic_regression(Xtemp,Ytemp); % Need to solve what is happening here if results are not satisfactory (order 6-7 particularly)
    catch
        temp = zeros(size(X,2)+1,1);
    end
    coeff = coeff + temp;
end
coeff = coeff./nSub;

end