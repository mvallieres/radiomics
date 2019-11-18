function [f,residuals] = corelucydenoiseNew(Y,H,DAMPAR22,wI,READOUT,SUBSMPL,idx,vec,num)
%CORELUCY Accelerated Damped Lucy-Richarson Operator.
%  Calculates function that when used with the scaled projected array 
%  produces the next iteration array that maximizes the likelihood that 
%  the entire suite satisfies the Poisson statistics. 
%
% See also DECONVLUCY and DECONVBLIND.

%  Copyright 1993-2003 The MathWorks, Inc.  
%  $Revision: 1.2.4.2 $ 

%   References
%   ----------
%   "Acceleration of iterative image restoration algorithms, by D.S.C. Biggs 
%   and M. Andrews, Applied Optics, Vol. 36, No. 8, 1997.
%   "Deconvolutions of Hubble Space Telescope Images and Spectra",
%   R.J. Hanisch, R.L. White, and R.L. Gilliland. in "Deconvolution of Images 
%   and Spectra", Ed. P.A. Jansson, 2nd ed., Academic Press, CA, 1997.

ReBlurred = real(ifftn(H.*fftn(Y)));

% 1. Resampling if needed
if SUBSMPL ~= 1,% Bin ReBlurred back to the sizeI for non-singleton dims
  
  %1.Reshape so that the-to-binned dimension separates into two
  %dimensions, with one of them consisting of elements of a single bin.
  ReBlurred = reshape(ReBlurred,vec);

  %2. Bin (==calculate mean) along the first of the-to-binned dimension,
  %that dimension consists of the bin elements. Reshape to get rid off
  for k = num,% new appeared singleton.
    vec(k) = [];
    ReBlurred = reshape(mean(ReBlurred,k),vec);
  end
  
end;


% Define Residual (This part we added, Andre Diamant)

Residual=wI-ReBlurred;
%fprintf('Residual calculated and mean equals %.5f \n',mean(Residual(:)))

% Apply the wavelet xform denoising method

Residual=denoise(Residual,3);
fprintf('DONE: Mean residual equals %.5f \n',mean(Residual(:)))

% 2. An Estimate for the next step
ReBlurred = ReBlurred + READOUT;
ReBlurred(ReBlurred == 0) = eps;
AnEstim = (ReBlurred+Residual)./ReBlurred + eps;


% 3. Damping if needed
if DAMPAR22 == 0,% No Damping
  ImRatio = AnEstim(idx{:});
else % Damping of the image relative to DAMPAR22 = (N*sigma)^2
  gm = 10;
  g = (wI.*log(AnEstim)+ ReBlurred - wI)./DAMPAR22;
  g = min(g,1);
  G = (g.^(gm-1)).*(gm-(gm-1)*g);
  ImRatio = 1 + G(idx{:}).*(AnEstim(idx{:}) - 1);
end;

f = fftn(ImRatio);
residuals = sum(Residual(:));
