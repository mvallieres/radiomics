function waitBatch(pathCheck,time,nBatch)
% -------------------------------------------------------------------------
% function waitBatch(pathCheck,time,nBatch)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function implements a waiting loop ensuring that all the
% computations from all parallel batch are done.
% -------------------------------------------------------------------------
% INPUTS:
% 1. pathCheck: Full path to the directory where the 'batch1_end',
%   'batch2_end', etc. (parallel checkpoints) are saved.
% 2. time: Number of seconds to wait before checking if parallel 
%          computations are done.
%          --> Ex: 60
% 3. nBatch: Number of parallel batch.
%            --> Ex: 8
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

startpath = pwd;

cd(pathCheck)
if nBatch
    while 1
        pause(time);
        check = zeros(nBatch,1);
        for i = 1:nBatch
            check(i) = exist(['batch',num2str(i),'_end'],'file');
        end
        if sum(check) == nBatch*2
            break
        end
    end
end

cd(startpath)
end