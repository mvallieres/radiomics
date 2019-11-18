function [patients] = batchPatients(nPatient,nBatch)
% -------------------------------------------------------------------------
% function [patients] = batchPatients(nPatient,nBatch)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function separates a given number of patients into a given number
% of batch for subsequent parallel computations.
% -------------------------------------------------------------------------
% INPUTS:
% - pathWORK: Full path to the HN WORKSPACE directory.
% - nPatient: Total number of patients.
% - nBatch: Number of parallel batch.
% -------------------------------------------------------------------------
% OUTPUTS
% - patients: Cell of 'nBatch' arrays, where each array specifies the
%             patients numbers in each batch.
% -------------------------------------------------------------------------
% AUTHOR(S): Martin Vallieres <mart.vallieres@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: July 2015
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


% FIND THE NUMBER OF PATIENTS IN EACH BATCH
patients = cell(1,nBatch);
if nBatch
    nP = nPatient / nBatch; nSup = ceil(nP); nInf = floor(nP);
    if nSup ~= nInf
        nSubInf = nBatch - 1; nSubSup = 1; total = nSubInf*nInf + nSubSup*nSup;
        while total ~= nPatient
            nSubInf = nSubInf - 1; nSubSup = nSubSup + 1;
            total = nSubInf*nInf + nSubSup*nSup;
        end
        nP = [repmat(nInf,[1,nSubInf]),repmat(nSup,[1,nSubSup])];
    else % The number of patients in all batches will be the same
        nP = repmat(nSup,[1,nBatch]);
    end
    start = 1;
    for i = 1:nBatch
        patients{i} = start:(start+nP(i)-1);
        start = start+nP(i);
    end
end

end