function [J,residuals] = deconvlucydenoiseNew(varargin) %CHANGED - added meanvec
%DECONVLUCY Deblur image using Lucy-Richardson method.
%   J = DECONVLUCY(I,PSF) deconvolves image I using Lucy-
%   Richardson algorithm, returning deblurred image J. The assumption is
%   that the image I was created by convolving a true image with a
%   point-spread function PSF and possibly by adding noise.
%   
%   I can be an N-Dimensional array.
%
%   To improve the restoration, additional parameters can be passed in
%   (use [] as a place holder if an intermediate parameter is unknown):
%   J = DECONVLUCY(I,PSF,NUMIT)
%   J = DECONVLUCY(I,PSF,NUMIT,DAMPAR)
%   J = DECONVLUCY(I,PSF,NUMIT,DAMPAR,WEIGHT)
%   J = DECONVLUCY(I,PSF,NUMIT,DAMPAR,WEIGHT,READOUT)
%   J = DECONVLUCY(I,PSF,NUMIT,DAMPAR,WEIGHT,READOUT,SUBSMPL), where
%
%   NUMIT   (optional) is the number of iterations (default is 10).
%
%   DAMPAR  (optional) is an array that specifies the threshold deviation
%   of the resulting image from the image I (in terms of the standard 
%   deviation of Poisson noise) below which the damping occurs. The 
%   iterations are suppressed for the pixels that deviate within the 
%   DAMPAR value from their original value. This suppresses the noise 
%   generation in such pixels, preserving necessary image details
%   elsewhere. Default is 0 (no damping).
%
%   WEIGHT  (optional) is assigned to each pixel to reflect its recording
%   quality in the camera. A bad pixel is excluded from the solution by
%   assigning it zero weight value. Instead of giving a weight of one for
%   good pixels, you can adjust their weight according to the amount of
%   flat-field correction. Default is a unit array of the same size as 
%   input image I.
%
%   READOUT (optional) is an array (or a value) corresponding to the
%   additive noise (e.g., background, foreground noise) and the variance 
%   of the read-out camera noise. READOUT has to be in the units of the
%   image. Default is 0.
%
%   SUBSMPL (optional) denotes subsampling and is used when the PSF is
%   given on a grid that is SUBSMPL times finer than the image. Default
%   is 1.
%
%   Note that the output image J could exhibit ringing introduced by the
%   discrete Fourier transform used in the algorithm. To reduce the
%   ringing use I = EDGETAPER(I,PSF) prior to calling DECONVLUCY.
%
%   Note also that DECONVLUCY allows you to resume deconvolution starting
%   from the results of an earlier DECONVLUCY run. To initiate this
%   syntax, the input image I has to be passed in as cell array, {I}.
%   Then the output J becomes a cell array and can be passed as the input
%   array into the next DECONVLUCY call. The input cell array can contain
%   one numeric array (on initial call), or four numeric arrays (when it
%   is the output from a previous run of DECONVLUCY). The output J
%   contains four elements, where J{1}=I, J{2} is the image resulted from
%   the last iteration, J{3} is the image from one before last iteration,
%   J{4} is an array used internally by the iterative algorithm.
%
%   Class Support
%   -------------
%   I and PSF can be uint8, uint16, int16, double, or single. DAMPAR and
%   READOUT must have the same class as the input image. Other inputs have to
%   be double. The output image (or the first array of the output cell) has
%   the same class as the input image.
%
%   Example
%   -------
%
%      I = checkerboard(8);
%      PSF = fspecial('gaussian',7,10);
%      V = .0001;
%      BlurredNoisy = imnoise(imfilter(I,PSF),'gaussian',0,V);
%      WT = zeros(size(I));WT(5:end-4,5:end-4) = 1;
%      J1 = deconvlucy(BlurredNoisy,PSF);
%      J2 = deconvlucy(BlurredNoisy,PSF,20,sqrt(V));
%      J3 = deconvlucy(BlurredNoisy,PSF,20,sqrt(V),WT);
%      subplot(221);imshow(BlurredNoisy);
%                     title('A = Blurred and Noisy');
%      subplot(222);imshow(J1);
%                     title('deconvlucy(A,PSF)');
%      subplot(223);imshow(J2);
%                     title('deconvlucy(A,PSF,NI,DP)');
%      subplot(224);imshow(J3);
%                     title('deconvlucy(A,PSF,NI,DP,WT)');
%
%   See also DECONVWNR, DECONVREG, DECONVBLIND, EDGETAPER, IMNOISE, PADARRAY, 
%            PSF2OTF, OTF2PSF.

%   Copyright 1993-2011 The MathWorks, Inc.
%   $Revision: 1.6.4.10 $
%

%   References
%   ----------
%   "Acceleration of iterative image restoration algorithms, by D.S.C. Biggs 
%   and M. Andrews, Applied Optics, Vol. 36, No. 8, 1997.
%   "Deconvolutions of Hubble Space Telescope Images and Spectra",
%   R.J. Hanisch, R.L. White, and R.L. Gilliland. in "Deconvolution of Images 
%   and Spectra", Ed. P.A. Jansson, 2nd ed., Academic Press, CA, 1997.

% Parse inputs to verify valid function calling syntaxes and arguments
[J,PSF,NUMIT,DAMPAR,READOUT,WEIGHT,SUBSMPL,sizeI,classI,numNSdim]=...
  parse_inputs(varargin{:});

% 1. Prepare PSF. If PSF is known at a higher sampling rate, it has to be
% padded with zeros up to sizeI(numNSdim)*SUBSMPL in all non-singleton
% dimensions. Or its OTF could take care of it:
sizeOTF = sizeI;
sizeOTF(numNSdim) = SUBSMPL*sizeI(numNSdim);
H = psf2otf(PSF,sizeOTF);

% 2. Prepare parameters for iterations
%
% Create indexes for image according to the sampling rate
idx = repmat({':'},[1 length(sizeI)]);
for k = numNSdim,% index replicates for non-singleton PSF sizes only
  idx{k} = reshape(repmat(1:sizeI(k),[SUBSMPL 1]),[SUBSMPL*sizeI(k) 1]);
end

wI = max(WEIGHT.*(READOUT + J{1}),0);% at this point  - positivity constraint
J{2} = J{2}(idx{:});
scale = real(ifftn(conj(H).*fftn(WEIGHT(idx{:})))) + sqrt(eps);
clear WEIGHT;
DAMPAR22 = (DAMPAR.^2)/2;

if SUBSMPL~=1,% prepare vector of dimensions to facilitate the reshaping
  % when the matrix is binned within the iterations.
  vec(2:2:2*length(sizeI)) = sizeI;
  vec(2*numNSdim-1) = -1;
  vec(vec==0) = [];
  num = fliplr(find(vec==-1));
  vec(num) = SUBSMPL;
else
  vec = [];    
  num = [];
end

% 3. L_R Iterations
% 
residuals=zeros(1,NUMIT); %CHANGED - initializing residualvector

lambda = 2*any(J{4}(:)~=0);
% h = waitbar(0,'Applying PVE correction to input volume ... please wait'); % waitbar
tic
for k = lambda + (1:NUMIT)
    
  % 3.a Make an image predictions for the next iteration    
  if k > 2,
    lambda = (J{4}(:,1).'*J{4}(:,2))/(J{4}(:,2).'*J{4}(:,2) +eps);
    lambda = max(min(lambda,1),0);% stability enforcement
  end
  Y = max(J{2} + lambda*(J{2} - J{3}),0);% plus positivity constraint
  
  % 3.b  Make core for the LR estimation
  fprintf('Iteration %u ... ',k)
  [CC,residuals(k)] = corelucydenoiseNew(Y,H,DAMPAR22,wI,READOUT,SUBSMPL,idx,vec,num); %CHANGED to "denoise" and to call meanvec
  
  % 3.c Determine next iteration image & apply positivity constraint
  J{3} = J{2};
  J{2} = max(Y.*real(ifftn(conj(H).*CC))./scale,0);  
  clear CC;
  J{4} = [J{2}(:)-Y(:) J{4}(:,1)];
  % waitbar(k/NUMIT,h,'Applying PVE correction to input volume ... please wait') % waitbar
end
% waitbar(1,h,'Applying PVE correction to input volume ... DONE'), pause(1) % waitbar
toc
% close(h)
clear wI H scale Y;

% 4. Convert the right array (for cell it is first array, for notcell it is
% second array) to the original image class & output whole thing
num = 1 + strcmp(classI{1},'notcell');
if ~strcmp(classI{2},'double'),
  J{num} = changeclass(classI{2},J{num});
end

if num==2,% the input & output is NOT a cell
  J = J{2};
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function: parse_inputs 
function [J,PSF,NUMIT,DAMPAR,READOUT,WEIGHT,SUBSMPL,sizeI,classI,numNSdim] = ...
    parse_inputs(varargin)
%
% Outputs:
% I=J{1}   the input array (could be any numeric class, 2D, 3D)
% PSF      operator that distorts the ideal image
% numNSdim non-singleton dimensions of PSF
%
% Defaults:
%
NUMIT = [];NUMIT_d = 10;% Number of  iterations, usually produces good
                        % result by 10.
DAMPAR =[];DAMPAR_d = 0;% No damping is default
WEIGHT =[];             % All pixels are of equal quality, flat-field is one
READOUT=[];READOUT_d= 0;% Zero readout noise or any other
           % back/fore/ground noise associated with CCD camera.
           % Or the Image is corrected already for this noise by user.
SUBSMPL= [];SUBSMPL_d= 1;% Image and PSF are given at equal resolution,
           % no over/under sampling at all.

narginchk(2,7);

% First, assign the inputs starting with the image
%
if iscell(varargin{1}),% input cell is used to resume interrupted iterations
  classI{1} = 'cell';% or interrupt the iteration to resume it later
  J = varargin{1};
else % no-cell array is used to do a single set of iterations
  classI{1} = 'notcell';  
  J{1} = varargin{1};% create a cell array in order to do the iterations
end;

% check the Image, which is the first array of the cell
classI{2} = class(J{1});

validateattributes(J{1},{'uint8' 'uint16' 'double' 'int16','single'},...
              {'real' 'nonempty' 'finite'},mfilename,'I',1);

if length(J{1})<2,
    error(message('images:deconvlucy:inputImagesMustHaveAtLeast2Elements'))
elseif ~isa(J{1},'double'),
    J{1} = im2double(J{1});
end

% now since the image is OK&double, we assign the rest of the J cell
len = length(J);
if len == 1,% J = {I} will be reassigned to J = {I,I,0,0}
  J{2} = J{1};
  J{3} = 0;
elseif len ~= 4,% J = {I,J,Jm1,gk} has to have 4 or 1 arrays
    error(message('images:deconvlucy:inputCellMustHave1or4Elements'));
else % check if J,Jm1,gk are double in the input cell
  if ~all([isa(J{2},'double'),isa(J{3},'double'),isa(J{4},'double')]),
    error(message('images:deconvlucy:inputImageCellElementsMustBeDouble'))
  end
end;

% Second, Assign the rest of the inputs:
%
PSF = varargin{2};%      deconvlucy(I,PSF)
switch nargin
case 3,%                 deconvlucy(I,PSF,NUMIT)
  NUMIT = varargin{3};
case 4,%                 deconvlucy(I,PSF,NUMIT,DAMPAR) CHANGED
  NUMIT = varargin{3};
  DAMPAR = varargin{4};
case 5,%                 deconvlucy(I,PSF,NUMIT,DAMPAR,WEIGHT)
  NUMIT = varargin{3};
  DAMPAR = varargin{4};
  WEIGHT = varargin{5};
case 6,%                 deconvlucy(I,PSF,NUMIT,DAMPAR,WEIGHT,READOUT)
  NUMIT = varargin{3};
  DAMPAR = varargin{4};
  WEIGHT = varargin{5};
  READOUT = varargin{6};
case 7,%                 deconvlucy(I,PSF,NUMIT,DAMPAR,WEIGHT,READOUT,SUBSMPL)
  NUMIT = varargin{3};
  DAMPAR = varargin{4};
  WEIGHT = varargin{5};
  READOUT = varargin{6};
  SUBSMPL = varargin{7};
end

% Third, Check validity of the input parameters: 
%
% NUMIT check number of iterations
if isempty(NUMIT),
  NUMIT = NUMIT_d;
else  
  validateattributes(NUMIT,{'double'},{'scalar' 'positive' 'finite'},...
                mfilename,'NUMIT',3);
end

% SUBSMPL check sub-sampling rate
if isempty(SUBSMPL),
  SUBSMPL = SUBSMPL_d;
else
  validateattributes(SUBSMPL,{'double'},{'scalar' 'positive' 'finite'},...
                mfilename,'SUBSMPL',7);
end

% PSF array
[sizeI, sizePSF] = padlength(size(J{1}), size(PSF));
numNSdim = find(sizePSF~=1);
if prod(sizePSF)<2,
  error(message('images:deconvlucy:psfMustHaveAtLeast2Elements'))
elseif all(PSF(:)==0),
  error(message('images:deconvlucy:psfMustNotBeZeroEverywhere'))
elseif any(sizePSF(numNSdim)/SUBSMPL > sizeI(numNSdim)),
  error(message('images:deconvlucy:psfMustBeSmallerThanImage'))
end
if length(J)==3,% assign the 4-th element of input cell now
  J{4}(prod(sizeI)*SUBSMPL^length(numNSdim),2) = 0;
end;

% DAMPAR check damping parameter
if isempty(DAMPAR),
  DAMPAR = DAMPAR_d;
elseif (numel(DAMPAR)~=1) && ~isequal(size(DAMPAR),sizeI),
  error(message('images:deconvlucy:damparMustBeSameSizeAsImage'))
elseif ~isa(DAMPAR,classI{2}),
  error(message('images:deconvlucy:damparMustBeSameClassAsInputImage'))
elseif ~strcmp(classI{2},'double'),
  DAMPAR = im2double(DAMPAR);
end

validateattributes(DAMPAR,{'double'},{'finite'},mfilename,'DAMPAR',4);

% READOUT check read-out noise
if isempty(READOUT),
  READOUT = READOUT_d;
elseif (numel(READOUT)~=1) && ~isequal(size(READOUT),sizeI),
  error(message('images:deconvlucy:readoutMustBeSameSizeAsImage'))
elseif ~isa(READOUT,classI{2}),
  error(message('images:deconvlucy:readoutMustBeSameClassAsInputImage'))
elseif ~strcmp(classI{2},'double'),
  READOUT = im2double(READOUT);
end

validateattributes(READOUT,{'double'},{'finite'},mfilename,'READOUT',6);

% WEIGHT check weighting
if isempty(WEIGHT),
  WEIGHT = ones(sizeI);
else
    validateattributes(WEIGHT,{'double'},{'finite'},mfilename,'WEIGHT',5);    
    if (numel(WEIGHT)~=1) && ~isequal(size(WEIGHT),sizeI),
      error(message('images:deconvlucy:weightMustBeSameSizeAsImage'))
    elseif numel(WEIGHT)== 1,
      WEIGHT = repmat(WEIGHT,sizeI);
    end
end
