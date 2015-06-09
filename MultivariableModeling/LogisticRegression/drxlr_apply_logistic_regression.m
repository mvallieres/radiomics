function [mu,b,PValuesColumn,WaldColumn,LR_chi2,convState]=drxlr_apply_logistic_regression(x,y,logiter_param,logtol_param)
%DREXLER subfunction 
%Written by Issam El Naqa 2003-2005
%Extracted for generalized use 2005, AJH
% logistic regression function
% LM: APA 07/13/2006, added convState as output parameter to represent the
% state of convergence
%
% Copyright 2010, Joseph O. Deasy, on behalf of the DREES development team.
% 
% This file is part of the Dose Response Explorer System (DREES).
% 
% DREES development has been led by:  Issam El Naqa, Aditya Apte, Gita Suneja, and Joseph O. Deasy.
% 
% DREES has been financially supported by the US National Institutes of Health under multiple grants.
% 
% DREES is distributed under the terms of the Lesser GNU Public License. 
% 
%     This version of DREES is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
% DREES is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
% See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with DREES.  If not, see <http://www.gnu.org/licenses/>.

if ~exist('logiter_param') | ~exist('logtol_param')  % if missing use default
    logiter_param=300;
    logtol_param=1e-6;
end
[mu, b,se, LR_chi2,convState]=drxlr_logistic_regression(x,y, logiter_param,logtol_param);
Wald=b(1:end-1)./se(1:end-1);
WaldColumn=Wald(:);
PValuesColumn=drxlr_get_p_gaussian(WaldColumn);
return

