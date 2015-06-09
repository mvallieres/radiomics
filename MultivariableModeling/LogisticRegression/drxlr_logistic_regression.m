function [mu, b, se, LR_chi2, convState] = drxlr_logistic_regression(x,y,iterlim,convcrit)
% regression by logistic model
% discrete reponse and different types of explanatory variables
% x: column vector of explanatory variables
% y: column vector of observations
% solved by maximum likelihood estimates using iterative  weighted
% least-squared.
% written by Issam El Naqa Spring 2003
% Extracted for generalized use 2005, AJH
% LM: APA 07/13/2006, added convState as output parameter to represent the
% state of convergence. convState is a 1x3 vector with 1st element being
% converged Yes or no, second element being number of iterations for
% convergence and third element being convergence criteria condition
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

% Initialize parameters
%warning('off');
x = [x ones(size(x,1),1)]; % constant factor
iter = 0;
[n p] = size(x);
seps = sqrt(eps);
eta = drxlr_rlogit((y + 0.5)/2); % initialize
b = zeros(p,1);
b0 = b+1;
% Start iterating
while(1)
    iter = iter+1;
    mu = drxlr_invlogit(eta);
    mu = max(0, min(1, mu));
    % Check stopping conditions
    if (~any(abs(b-b0) > convcrit * max(seps, abs(b0))))
        convState = [1 iter-1 max(abs(b-b0))];
        break;
    end
    if (iter>iterlim)
        warning('drxlr_logistic_regression: Number of iterations exceeded without convergence (convergence tolerance = %e', convcrit)
        convState = [2 iter-1 max(abs(b-b0))];
        break;  % iteration limit reached, exit!
    end
    deta = drxlr_d_logit(mu);
    z = eta + (y - mu) .* deta;
    vy=drxlr_binomial_variance(mu);
    w = 1 ./ max(eps, (deta .^ 2) .* vy);
    b0 = b;
    [b,R,se] = drxlr_wfit(z, x, w, p);
    eta = x * b;
end

%%% calculate likelihood ratio test
LR_chi2=drxlr_get_lrt(y,mu);
return
