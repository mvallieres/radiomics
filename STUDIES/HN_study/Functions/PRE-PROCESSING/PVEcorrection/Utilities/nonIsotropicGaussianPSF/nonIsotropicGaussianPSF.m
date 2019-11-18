function [outKer] = nonIsotropicGaussianPSF(inSigma,varargin)
%% NONISOTROPICGAUSSIANPSF - Creates a isotropic / non isotropic 3D Gaussian kernel.
%Uses the Statistical toolbox if installed (mvnpdf). Works without Statistical Toolbox, but slower and
% with higher limitation when it comes to the size of the PSF support..
%
% Syntax:  [outKer]     = nonIsotropicPSF(inSigma)
%                                  = nonIsotropicPSF(inSigma,sizeDomain)
%                                  = nonIsotropicPSF(inSigma,sizeDomain,precision)
% Inputs:
%    inSigma - scalar (isotropic) or 3x1 vector with the standard deviation of the Gaussian kernel (reminder :
%    sigma=FWHM/(2*sqrt(2*ln(2))) )
%   sizeDomain (optional) - Define the size of the support of the PSF (default: 2.1*max(inSigma))
%   precision (optional) - Add 'single' as input to calculate in single precision (default: double).
%   The precision is important when using the function without the Statistical Toolbox.
%
% Outputs:
%    outKer - 3D non isotropic Gaussian kernel
%
% Example:
%   outKer = nonIsotropicGaussianPSF([5.12 5.9 5.8],3);
%   If Out of Memory, try : 
%   outKer = nonIsotropicGaussianPSF([5.12 5.9 5.8],3,'single');
%
% Other m-files required: areTheseToolboxesInstalled (added in this file)
% Subfunctions: none
% MAT-files required: none

% Author: Christopher Coello
% Work address: Preclinical PET/CT Unit
% Email address: s.c.coello@medisin.uio.no
% Website: http://www.med.uio.no/imb/english/services/public/pet/
% 2012/03/13; Last revision: 2012/03/16
% Created with Matlab version: 7.13.0.564 (R2011b)
%%
% Init
outKer=0;
% Input check
switch length(varargin),
    case 2
        if strcmp(varargin{2},'single'),
            indS=1;
        else
            display('Wrong precision argument. Only single is accepted (default : double)');return;
        end
        extentSupport=varargin{1};
    case 1
        if isnumeric(varargin{1})
            extentSupport=varargin{1};
            indS=0;
        elseif strcmp(varargin{1},'single'),
            indS=1;
            extentSupport=2.1;
        else
            display('Wrong precision argument. Only single is accepted (default : double)');return;
        end
    otherwise
        extentSupport=2.1;
        indS=0;
end

% Defines the covariance matrix
if isscalar(inSigma),
    inSigma=ones(3,1)*inSigma;
end
if indS,
    S=single(diag(inSigma));
else
    S=diag(inSigma);
end

% Defines the size of the kernel in relation with the Gaussian properties
supportSize=ceil(extentSupport*max(inSigma));

% Initialise the volume where the kernel is to be drawn
[X,Y,Z]=meshgrid(-supportSize:1:supportSize);

% Put it in columns
if indS,
    coord=single([X(:) Y(:) Z(:)]);
else
    coord=[X(:) Y(:) Z(:)];
end
mu=mean(coord);

% Define the output size
outsize=size(X);


if (areTheseToolboxesInstalled('Statistics Toolbox')),
    
    %with statistics toolbox
    clear X Y Z
    
    % Check mvnpdf help for details on this function
    p = mvnpdf(coord, mu, S);
    
    % Output vector is reshaped to volume.
    outKer=reshape(p,outsize);
    
else
    
    %without statistics toolbox
    
    % Calculates the amplitude of the PSF
    if indS,
        AmplitudePSF=single(1/((2*pi)^(3/2)*det(S)^(1/2)));
    else
        AmplitudePSF=1/((2*pi)^(3/2)*det(S)^(1/2));
    end
    
    % Center the volume
    clear coord
    if indS,
        c_coord=single([X(:)-mu(1) Y(:)-mu(2) Z(:)-mu(3)]);
    else
        c_coord=[X(:)-mu(1) Y(:)-mu(2) Z(:)-mu(3)];
    end
    
    % Apply the equation of the Gaussian distribution
    T_hand=AmplitudePSF*exp(-0.5*c_coord*(S^-1)*c_coord');
    
    % Reshape to volume
    outKer=reshape(diag(T_hand),outsize);
    
end

% Check if the support is big enough
if ((1-sum(outKer(:)))*100>0.01),
    display(['Warning : ' num2str((1-sum(outKer(:)))*100) '% of the Gaussian is outside the domain.'])
    display('If not acceptable, widen the domain of the PSF by defining a superior second argument (default: 2.1*max(inSigma)).')
end


function tf = areTheseToolboxesInstalled(requiredToolboxes)
%ARETHESETOOLBOXESINSTALLED takes a cell array of toolbox names and checks whether they are currently installed
% SYNOPSIS tf = areTheseToolboxesInstalled(requiredToolboxes)
%
% INPUT requiredToolboxes: cell array with toolbox names to test for. Eg.
%        {'MATLAB','Image Processing Toolbox'}
%
% OUTPUT tf: true or false if the required toolboxes are installed or not
%%%%%%%%%%%%%%%%%%%%%%%%%%

% get all installed toolbox names
v = ver;

% collect the names in a cell array
[installedToolboxes{1:length(v)}] = deal(v.Name);

% check
tf = all(ismember(requiredToolboxes,installedToolboxes));