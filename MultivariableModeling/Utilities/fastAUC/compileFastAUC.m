function compileFastAUC(OS)
% -------------------------------------------------------------------------
% function compileFastAUC(OS)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function finds the full path to the fastAUC.cpp file, and then
% compiles it. The code fastAUC.cpp is written by Enric Junqué de Fortuny
% and is available at: <http://www.mathworks.com/matlabcentral/fileexchange/41258-faster-roc-auc>
% -------------------------------------------------------------------------
% INPUTS:
% - OS: String specifying the type of operating system. 
%       --> Now supports 'Linux'
%       --> Future release: 'Unix', 'Windows'
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

startpath = pwd;

% Find fastAUC.cpp in the system
if strcmp(OS,'Linux')
   [~,~] = system('find / -name "fastAUC.cpp">temp.txt');
end

% Reads in a temporary file where fastAUC.cpp is located
fid = fopen('temp.txt');
temp = fgetl(fid);
pathFastAUC = temp(1:end-12);

% Clean up temporary file
fclose(fid);
delete('temp.txt')

% Go to fastAUC.cpp directory, and start compilation
cd(pathFastAUC)
mex fastAUC.cpp
addpath(pathFastAUC), savepath

cd(startpath)
end
