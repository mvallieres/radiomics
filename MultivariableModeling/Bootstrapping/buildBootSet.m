function [bootSam,testSets] = buildBootSet(Y,nBoot,IABR)
% -------------------------------------------------------------------------
% function [bootSam] = buildBootSet(outcome,nBoot,IABR)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function creates a matrix of 'nBoot' bootstrap samples. The matrix
% is composed of instance indexes for each bootstrap sample. In addition,
% if specified by the user, the function implements the 'imbalance-adjusted 
% bootstrap resampling'method: it duplicates the number of positive instance
% by a factor equal to the number of negative instances, and it duplicates
% the number of negative instance by a factor equal to the number of 
% positive instances, such that the probability of picking a positive and a 
% negative instance in the bootstrap sample is equal. See ref. [1] for more
% details.
% -------------------------------------------------------------------------
% REFERENCE:
% [1] Vallieres, M. et al. (2015). A radiomics model from joint FDG-PET and 
%     MRI texture features for the prediction of lung metastases in soft-tissue 
%     sarcomas of the extremities. Physics in Medicine and Biology, 60(14), 
%     5471-5496. doi:10.1088/0031-9155/60/14/5471
% -------------------------------------------------------------------------
% INPUTS:                             
% - Y: Column vector of size [nInst X 1] specifying the outcome status 
%      (1 or 0) for all instances.
% - nBoot: Number of bootstrap samples.
% - IABR: (optional argument). Use 'adjust' to perform imbalance-adjusted
%         bootstrap resampling (IABR), and use no argument for regular 
%         resampling.
% -------------------------------------------------------------------------
% OUTPUTS:
% - bootSam: Matrix of size [nInst X nBoot], where 'nInst' is the number of
%            instance in Y.
% - testSets: Cell of size [1 X nBoot], where each entry is a testing set
%             vector (Testing sets for different bootstrap samples may not 
%             be of the same length).
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
end

indPos = find(Y); % indexes of postivie instances
indNeg = find(~Y); % indexes of negative instances
nNeg = length(indNeg); % Number of negative instances
nPos = length(indPos); % Number of positive instances
nInst = length(indPos) + length(indNeg); % Total number of instances

if nargin == 3 && strcmp(IABR,'adjust') % Adjust probability of picking positive and negative instances
    indPos = repmat(indPos,nNeg,1);   
    indNeg = repmat(indNeg,nPos,1);
end
ind = [indPos;indNeg];
ind = shufflePartition(ind);
ind = shufflePartition(ind);

% Obtaining the bootstrap samples
nSize = length(ind);
bootSam = ind(ceil(nSize .* rand(nInst,nBoot)));
testSets = findBootTestSet(bootSam);

% Verification for each bootstrap sample: Both training and testing sets 
% must have at least 2 instances of each class, and these instances must be 
% different.
for n = 1:nBoot
    training = bootSam(:,n);
    testing = testSets{n};
    while ( sum(Y(testing))<2 || sum(Y(testing))>(length(Y(testing))-2) || sum(Y(training))<2 || sum(Y(training))>(length(Y(training))-2) || length(unique(training))==1 || length(unique(testing))==1 )
        bootSam(:,n) = ind(ceil(nSize .* rand(nInst,1)));
        testSets(n) = findBootTestSet(bootSam(:,n));
        training = bootSam(:,n);
        testing = testSets{n};
    end
end

end
