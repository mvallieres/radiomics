function [pathMINE] = findMINE(OS)
% -------------------------------------------------------------------------
% function [pathMINE] = findMINE(OS)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function finds the full path to the MINE.jar executable. The
% executable can be downloaded at: <http://www.exploredata.net/>.
% -------------------------------------------------------------------------
% INPUTS:
% - OS: String specifying the type of operating system. 
%       --> Now supports 'Linux'
%       --> Future release: 'Unix', 'Windows'
% -------------------------------------------------------------------------
% OUTPUTS:
% -  pathMINE: Full path to the MINE.jar executable.
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


% Find MINE.jar in the system
if strcmp(OS,'Linux')
   [~,~] = system('find / -name "MINE.jar">temp.txt');
end

% Reads in a temporary file where MINE.jar is located
fid = fopen('temp.txt');
temp = fgetl(fid);
pathMINE = temp(1:end-9);

% Clean up temporary file
fclose(fid);
delete('temp.txt')

end