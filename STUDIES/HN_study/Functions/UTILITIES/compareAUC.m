function [p,CI] = compareAUC(X1,X2,Y)
% -------------------------------------------------------------------------
% function [p,CI] = compareAUC(X1,X2,Y)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function calculates the significance (p-value) for the difference
% between two CORRELATED AUCs (i.e. computed from the same samples) 
% computed for two different prediction situations. The implementation of 
% this function is based on ref. [1].
% -------------------------------------------------------------------------
% REFERENCE:
% [1] DeLong, E. R., DeLong, D. M. & Clarke-Pearson, D. L. Comparing the 
%     areas under two or more correlated receiver operating characteristic 
%     curves: a nonparametric approach. Biometrics 44, 837–845 (1988).
% -------------------------------------------------------------------------
% INPUTS:
% 1. X1: [N X 1] vector of probabilities of outcome Y for situation 1.
% 2. X2: [N X 1] vector of probabilities of outcome Y for situation 2.
% 3. Y:  [N X 1] vector defining if each nth intance has the outcome or not
%        (1's or 0's), where N is the total number of patients.
% -------------------------------------------------------------------------
% OUTPUTS:
% 1. p: p-value for the difference in AUC.
% 2. CI: 95 % confidence interval for the difference in AUCs (diff), to be 
%    subsequently applied as diff - CI and diff + CI.
% -------------------------------------------------------------------------
% AUTHOR(S): Martin Vallieres <mart.vallieres@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: January 2017
% -------------------------------------------------------------------------
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


% SANITY CHECKS
if numel(X1) ~= numel(Y) || numel(X2) ~= numel(Y)
    error('The number of instances in vectors X1, X2, and Y must be the same')
end
Y = logical(Y);


% PARAMETERS OR INTEREST
sigma = 1.96; % For 95 % confidence interval
L = [1,-1]; % Contrast vector. We are always comparing 2 variables at a time, not more. This is hardcoded for now and should stay as is.


% INITIALIZATION
X = [X1,X2]; % Combining the variables into a single matrix. We are always comparing 2 variables at a time, not more.
nVar = size(X,2); % Number of variables we compare (referred to "r" in ref. [1], so here, "v" --> "r" in ref. [1])
nInst = numel(Y); nPos = sum(Y); nNeg = nInst - nPos; % Total number of instances, number of positive instances ("nPos" --> "m" in ref. [1]) and number of negative instances ("nNeg" --> "n" in ref. [1]), respectively
Vpos = zeros(nPos,nVar); Vneg = zeros(nNeg,nVar); % Matrices "Vpos" --> "V10" and "Vneg" --> "V01" in ref. [1], respectively. 
auc = zeros(1,nVar); % AUC estimates ("auc" --> "theta" in ref. [1])


% COMPUTATION OF AUC ESTIMATES ("auc" --> "theta"), VPOS ("Vpos" --> "V10") and VNEG ("Vneg" --> "V01")
Xpos = X(Y,:); Xneg = X(~Y,:); % Separating positive and negative instances ("Xpos" --> "X" and Xneg --> "Y" in ref. [1]. Here, our "Y" refers to the outcomes, or targets, to defined if an instance as the condition or not)
for v = 1:nVar
    for i = 1:nPos
        val = Xpos(i,v);
        Vpos(i,v) = sum(Xneg(:,v) < val) + 0.5*sum(Xneg(:,v) == val);
    end
    for i = 1:nNeg
        val = Xneg(i,v);
        Vneg(i,v) = sum(Xpos(:,v) > val) + 0.5*sum(Xpos(:,v) == val);
    end
end
auc = sum(Vpos)/(nNeg*nPos);
Vpos = Vpos/nNeg; Vneg = Vneg/nPos;


% COMPUTATION OF COVARIANCE MATRICES ("Spos" --> "S10" and "Sneg" --> "S01" in ref. [1])
Spos = (1/nPos) * ((Vpos'*Vpos) - nPos*(auc*auc'));
Sneg = (1/nNeg) * ((Vneg'*Vneg) - nNeg*(auc*auc'));
S = (1/nPos)*Spos + (1/nNeg)*Sneg; 


% COMPUTATION OF CHI-SQUARE STATISTICS
LSL = L*S*L';
val_chi2 = (auc*L')*(inv(LSL))*(L*auc')*sqrt(2);
df_chi2 = rank(LSL);
p = 1-chi2cdf(val_chi2,df_chi2);
CI = sigma * sqrt(LSL);

end