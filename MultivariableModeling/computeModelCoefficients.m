function [coeff,resp,modelCI,coeff_boot] = computeModelCoefficients(X,Y,imbalance,seed)
% -------------------------------------------------------------------------
% function [coeff,resp,modelCI] = computeModelCoefficients(X,Y,imbalance,seed)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function computes the final model logistic regression coefficients 
% using bootstrap resampling, and it associated multivariable model 
% response and bootstrap confidence intervals. See ref. [1] for more 
% details. This function uses logistic regression utilities from DREES 
% <http://www.cerr.info/drees>, as well as a rounding function written by 
% Francois Beauducel available at: 
% <http://www.mathworks.com/matlabcentral/fileexchange/26212-round-with-significant-digits>
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
% - imbalance: String specifying the type of imbalance-adjustement strategy
%              employed. Either 'IABR' for imbalance-adjusted bootstrap
%              resampling (see ref.[1]), or 'IALR' for imbalance-adjusted
%              logistic regression (formal reference to come).
% - seed: (optional input). Numerical number to use as random generator 
%         seed for bootstrapping experiments.
%         --> Ex: 54288
% -------------------------------------------------------------------------
% OUTPUTS:
% - coeff: Column vector of size [nCoeff+1 X 1] specifying the final 
%          logistic regression coefficients. Last entry specify the offset
%          of the model.
% - resp: Column vector of size [nInst X 1] specifying the linear 
%         multivariable model response when 'coeff' is applied to 'X'.
% - modelCI: Column vector of size [nInst X 2] specifying the 95%
%            confidence interval on the multivariable model response as
%            defined by the 2.5 (modelCI(i,1)) and the 97.5 (modelCI(i,2))
%            percentiles, for the ith instance. See ref. [1] for more 
%            details.
% - coeff_boot: Array of size [nCoeff X 1000], where 1000 is the number of 
%               bootstrap samples used to calculate the mean coefficients 
%               ('coeff' output #1). Each column thus specify the logisitic 
%               regression model coefficients for a particular bootstrap
%               training sample. These set of bootstrap coeffcients could
%               therafter be used to calculate the confidence intervals on
%               predictions in new testing data.
% -------------------------------------------------------------------------
% AUTHOR(S): 
% - Martin Vallieres <mart.vallieres@gmail.com>
% - DREES development team <http://www.cerr.info/drees> (logistic regression)
% - Francois Beauducel (roundsd.m)
% -------------------------------------------------------------------------
% HISTORY:
% - Creation - May 2015
% - Revision I - July 2015: (including imbalance-adjusted logistic regression) 
% - Revision II - December 2016: Including initial seed as input (optional)
%                                for reproducibility of bootstrapping experiments
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
% --> Copyright (c) 2015, François Beauducel
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
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN 
% -------------------------------------------------------------------------


% INITIALIZATION
warning off
nBoot = 1000;
alpha = 0.05;
bound = 1; % One standard error
order = size(X,2);
coeff = zeros(order+1,nBoot);
tol = 10^6; % To avoid extremely large coefficients
if strcmp(imbalance,'IABR')
    logisticRegression = @(x,y) applyStandardLR(x,y);
    adjust = 'adjust';
elseif strcmp(imbalance,'IALR')
    logisticRegression = @(x,y) applyEnsembleLR(x,y);
    adjust = 'NoAdjust';
end


% RANDOM NUMBER GENERATOR SEED
if nargin == 4
    rng(seed)
else
    if ~RandStream.getGlobalStream.Seed
        rng('shuffle')
    end
end



% COMPUTING OVER ALL BOOTSTRAP SAMPLES
respBoot = zeros(size(X,1),nBoot);
modelCI = zeros(size(X,1),2);
for n = 1:nBoot
    average = inf;
    while average > tol
        [bootSam,~] = buildBootSet(Y,1,adjust);
        Xtrain = X(bootSam,:); Ytrain = Y(bootSam,1);
        coeff_temp = logisticRegression(Xtrain,Ytrain);
        coeff_temp(isnan(coeff_temp)) = 0;
        coeff(:,n) = coeff_temp;
        average = mean(abs(coeff(:,n)));
    end
    [respBoot(:,n)] = responseLR(X,coeff(:,n));
end
SE_coeff = bound.*(std(coeff')')./sqrt(nBoot); % Not used at the moment
coeff_boot = coeff;
coeff = mean(coeff')';


% CALCULATING THE BOOTSTRAP CONFIDENCE INTERVALS OF THE FINAL MODEL
for i = 1:size(X,1)
    modelCI(i,1)= prctile(respBoot(i,:),alpha/2*100);
    modelCI(i,2)= prctile(respBoot(i,:),(1-alpha/2)*100);
end


% ROUNDING COEFFICIENTS

% According to their standard errors, as in ref. [1]
% Note: This type of rounding seems too strong and may not be desirable, as it significantly affects the specificity of models
%       --> kept here as comments before further investigations
% for i = 1:order+1
%     if coeff(i) > 0
%         coeff(i) = roundsd(coeff(i),ceil(log10(abs(coeff(i)/roundsd(SE_coeff(i),1)))) + 1);
%     else
%         coeff(i) = roundsd(coeff(i),ceil(log10(abs(coeff(i)/roundsd(SE_coeff(i),1)))));
%     end
% end

% Rounding to 4 significant digits
% --> seems to preserve the predictive properties of models
coeff = roundsd(coeff,4);


% MULTIVARIABLE RESPONSE OF THE FINAL MODEL
[resp] = responseLR(X,coeff);

end