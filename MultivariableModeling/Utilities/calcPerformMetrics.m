function [AUC,sensitivity,specificity,accuracy] = calcPerformMetrics(response,Y,thresh)
% -------------------------------------------------------------------------
% function [AUC,sensitivity,specificity,accuracy] = calcPerformMetrics(response,outcome,thresh)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function computes AUC, sensitivity, specificity and accuracy metrics
% given a modle response, outcome and threshold. This function uses a fast 
% implementation of AUC calculation by Enric Junqué de Fortuny that is
% available at: <http://www.mathworks.com/matlabcentral/fileexchange/41258-faster-roc-auc>
% -------------------------------------------------------------------------
% INPUTS:
% - response: Column vector of size [nInst X 1], specifying the 
%             multivariable response for all instances (nInst is the total 
%             number of instances).
% - Y: Column vector of size [nInst X 1] specifying the outcome status 
%      (1 or 0) for all instances.
% - thresh: Threshold used to calculate sensitivity, specificity and
%           accuracy metrics (typically 0 for logistic regression).
% -------------------------------------------------------------------------
% OUTPUTS:
% - AUC
% - Sensitivity
% - Specificity
% - Accuracy
% -------------------------------------------------------------------------
% AUTHOR(S): 
% - Martin Vallieres <mart.vallieres@gmail.com>
% - Enric Junqué de Fortuny (fastAUC.cpp)
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


% Computing AUC
try
    AUC = fastAUC(Y,response,1);
catch
    try
        compileFastAUC('Linux')
        AUC = fastAUC(Y,response,1);
    catch
        try
            [~,~,~,AUC] = perfcurve(Y,response,1);
        catch
            AUC = 0.5;
            fprintf('\nSomething went wrong with the AUC calculations\n')
        end 
    end
end

% Classifications
TP = sum(response>=thresh & Y==1);
TN = sum(response<thresh & Y==0);
FP = sum(response>=thresh & Y==0);
FN = sum(response<thresh & Y==1);

% Computing performance metrics
sensitivity = TP/(TP + FN);
specificity = TN/(TN + FP);
accuracy = (TP + TN)/(TP + TN + FP + FN);

% Validation check of metrics
AUC(isnan(AUC)) = 0.5;
sensitivity(isnan(sensitivity)) = 0;
specificity(isnan(specificity)) = 0;
accuracy(isnan(accuracy)) = 0;

end