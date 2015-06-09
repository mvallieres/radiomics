function [valMIC] = applyMIC(pathMINE,variable,matData,jobName)
% -------------------------------------------------------------------------
% function [valMIC] = applyMIC(pathMINE,variable,matData,jobName)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function computes the Maximal Information Coefficient (MIC) between
% the feature contained in the vector 'variable' and all the features
% contained in the matrix 'matData'. MIC is computed using the executable
% MINE.jar, which can be downloaded at: <http://www.exploredata.net/>.
% -------------------------------------------------------------------------
% INPUTS:
% - pathMINE: Full path to the MINE.jar executable.
% - variable: Column cector of size [nInst X 1], where 'nInst' refers to the 
%             number of instances for the given feature tested.
% - matData: Matrix of size [nInst X nFeat], where 'nFeat' is the number of
%            features to be tested against 'variable'.
% - jobName: String specifying the name of the job to be sent to MINE.jar.
% -------------------------------------------------------------------------
% OUTPUTS:
% -  valMIC: Column vector of size [nFeat X 1], specifying the MIC between
%            each feature in 'matData' and the 'variable' features. Entry
%            number in 'valMIC' corresponds to the column number in
%            'matData'.
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
%
%    _______________________________________________________________
%
% MINE version 1.0.1d
% Copyright 2011 by David Reshef and Yakir Reshef.
%
% This application is licensed under a Creative Commons
% Attribution-NonCommercial-NoDerivs 3.0 Unported License.
%
% See http://creativecommons.org/licenses/by-nc-nd/3.0/ for more information.
% -------------------------------------------------------------------------


startpath=pwd;
cd(pathMINE)

% Writing the MIC input file
if isempty(strfind(jobName,'.csv'))
    jobName = [jobName,'.csv'];
end
id = (0:size(matData,2))';
data = [variable,matData]';
csvwrite(jobName,[id,data],0,0)

% Executing MIC
commandSys = ['java -jar MINE.jar ',jobName,' -masterVariable 0'];
[~,~] = system(commandSys);

% Reading the output from MIC
commandRead = [jobName,',mv=0,cv=0.0,B=n^0.6,Results.csv'];
micOutput = csvread(commandRead,1,1,[1 1 size(matData,2) 2]);

% Sorting the results
[val,ind] = sort(micOutput(:,1));
valMIC = val;
valMIC(1:end) = micOutput(ind(1:end),2);

% Cleaning up
delete(jobName)
delete(commandRead)
commandLast=[jobName,',mv=0,cv=0.0,B=n^0.6,Status.txt'];
delete(commandLast)

cd(startpath)
end