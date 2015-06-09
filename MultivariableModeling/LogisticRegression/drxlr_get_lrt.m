function LR_chi2=drxlr_get_lrt(y,mu)
% compute the likelihood ratio test 
% written by Issam El Naqa
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

seps = sqrt(eps);
n=length(y);
n1=sum(y); n0=sum(1-y);
Ls=n1*log(n1)+n0*log(n0)-n*log(n); % saturated
L=sum(y.*log(mu+seps)+(1-y).*log(1-mu+seps)); % fitted
LR_chi2=2*(L-Ls); % follows centered chi2 stats
return