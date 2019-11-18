function [significance] = benjamini_hochberg(pValues,FDR)
% - pvalues: Must be a vector of numerical values
% - FDR: False Discovery rate, between 0 and 1
% - OUTPUT: Vector of the same size as pValues specifying if each value is
% significant (1) or not (0).


% Consistency check
if ~isnumeric(FDR)
    error('FDR must be numeric')
else
    if FDR < 0 || FDR > 1
        error('FDR must be between 0 and 1')
    end
end

% Initialization
nPval = numel(pValues);
significance = zeros(nPval,1);

% Sorting the pvalues
[pSort,indSort] = sort(pValues);

% Critical value
rank = (1:nPval)';
cValue = rank/nPval*FDR;

% Finding the first significant p-value
valid = pSort < cValue;
indLast = find(valid,1,'last');

% Recording significant results
if ~isempty(indLast)
    valid(1:indLast) = 1;
    significance(indSort) = valid(1:end);
end

end