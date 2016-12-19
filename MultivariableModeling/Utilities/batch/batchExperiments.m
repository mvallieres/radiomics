function [param] = batchExperiments(setNames,outcomes,nBatch)
% -------------------------------------------------------------------------
% function [param] = batchExperiments(setNames,outcomes,nBatch)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function separates a given number of experiments into a given number
% of batch for subsequent parallel computations.
% -------------------------------------------------------------------------
% INPUTS:
% 1. setNames: Cell of strings specifying the name of feature sets.
%              --> Ex: {'set1','set2'}
% 2. outcomes: Structure of vectors specifying the status (1 or 0) for 
%              different cancer outcomes. Contains: outcomes.(nameOutcome1), 
%              outcomes.(nameOutcome2), where each outcome entry is a vector
%              of [nPatient X 1].
% 3. nBatch: Number of parallel batch.
%            --> Ex: 8
% -------------------------------------------------------------------------
% OUTPUTS
% 1. param: Cell vector of length equal to the number of outcomes. Each cell
%           contains a cell specifying the parameter numbers for each
%           experiment in each batch. The rows corresponds to the different
%           experiments, the first column to the feature set type, and the
%           second column the outcome type.
% -------------------------------------------------------------------------
% AUTHOR(S): Martin Vallieres <mart.vallieres@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: December 2016
%--------------------------------------------------------------------------
% STATEMENT:
% This file is part of <https://github.com/mvallieres/radiomics/>, 
% a package providing MATLAB programming tools for radiomics analysis.
% --> Copyright (C) 2015-2016  Martin Vallieres
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

nameOutcomes = fieldnames(outcomes);
nOutcomes = length(nameOutcomes);
nPatient = length(outcomes.(nameOutcomes{1}));
outWeight = zeros(1,nOutcomes);
for i = 1:nOutcomes
    temp = outcomes.(nameOutcomes{i});
    outWeight(i) = nPatient - sum(temp(~isnan(temp)));
end
outWeight = abs(outWeight);

type = 1:numel(setNames); nType = numel(type);
total = nType*nOutcomes;
param = cell(1,nBatch);
count = 0;
for i = 1:nBatch
    param{i} = [];
end
while count ~= total
    for i = 1:nBatch % From the start
        [~,ind] = max(outWeight);
        if isempty(type)
            type = 1:nType;
            outWeight(ind) = 0;
            [~,ind] = max(outWeight);
        end
        param{i} = [param{i};type(1),ind];
        count = count + 1;
        if count == total
            break
        end
        type(1) = [];
    end
    if count ~= total
        for i = 1:nBatch % From the end
            [~,ind] = max(outWeight);
            if isempty(type)
                type = 1:nType;
                outWeight(ind) = 0;
                [~,ind] = max(outWeight);
            end
            param{nBatch-i+1} = [param{nBatch-i+1};type(1),ind];
            count = count + 1;
            if count == total
                break
            end
            type(1) = [];
        end
    else
        break
    end
end

for i = 1:nBatch
    if isempty(param{i})
        param(i:end) = [];
        break
    end
end

temp = param; param = cell(1,numel(param));
for i = 1: numel(temp)
    nCell = size(temp{i},1);
    param{i} = cell(nCell,1);
    for j = 1:nCell
        param{i}{j} = {setNames{temp{i}(j,1)},nameOutcomes{temp{i}(j,2)}};
    end
end

end