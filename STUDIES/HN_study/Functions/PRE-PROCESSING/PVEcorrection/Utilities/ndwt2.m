function varargout = ndwt2(X,level,varargin)
%NDWT2 Non-decimated 2-D wavelet transform.
%   NDWT2 performs a multilevel 2-D non-decimated wavelet 
%   decomposition with respect to either a particular wavelet 
%   ('wname', see WFILTERS for more information) or particular  
%   wavelet filters you specify, and using a specified  
%   DWT extension mode (see DWTMODE).
%
%   WT = NDWT2(X,N,'wname','mode','ExtM') returns a structure 
%   which contains the non-decimated wavelet transform of the 
%   matrix X at the level N, N must be a strictly positive integer
%   (see WMAXLEV), 'wname' is a string containing the wavelet 
%   name and 'ExtM' is a string containing the extension mode.
%   WT = NDWT2(X,N,'wname') uses the default extension mode: 'sym'.
%
%   WT is a structure with the following fields:
%     sizeINI: contains the size of the 2-D array X.
%     level:   contains the level of the decomposition.
%     mode:    contains the name of the wavelet transform extension mode.
%     filters: is a structure with 4 fields LoD, HiD, LoR, HiR which
%              contain the filters used for DWT.
%         dec: is a 1 by (3*level+1) cell array containing the coefficients 
%              of the decomposition. dec{1} contains the coefficients of
%              the approximation and dec{j} (j = 2 to 3*level+1), contains  
%              the coefficients of the details from level level to the 
%              level 1, three details by level (LH, HL and HH where L
%              stands for low and H for high).
%       sizes: is a (level+1) by 2 array containing the size of the 
%              components. 
%
%   Instead of a single wavelet, you may specify two wavelets (i.e. one
%   wavelet per direction):
%   WT = NDWT2(X,W,...) with W = {'wname1','wname2'} or W a 
%   structure with 2 fields 'w1', 'w2' containing strings which
%   are the names of wavelets.
%
%   Instead of wavelets you may specify filters: 4 filters (2 for
%   decomposition and 2 for reconstruction) or 2x4 filters (one 
%   quadruplet per direction): WT = NDWT2(X,WF,...),
%   where WF must be a cell array (1x4) or (2x4): {LoD,HiD,LoR,HiR},
%   or a structure with the four fields 'LoD', 'HiD', 'LoR', 'HiR'.
%
%   Examples:
%       load noiswom; 
%       W1 = ndwt2(X,2,'db1')
%       W2 = ndwt2(X,3,'db1','mode','per')
%       W3 = ndwt2(X,3,{'db1','db2'},'mode','sym')
%       WF = W3.filters
%       W4 = ndwt2(X,3,WF,'mode','sym')
%
%   See also DWTMODE, INDWT2, WAVEINFO.

%   M. Misiti, Y. Misiti, G. Oppenheim, J.M. Poggi 09-Dec-2008.
%   Last Revision: 06-Feb-2011.
%   Copyright 1995-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.5 $ $Date: 2011/04/12 20:43:34 $

% Check arguments.
nbIn = length(varargin);
error(nargchk(1,5,nbIn,'struct'))

LoD = cell(1,2); HiD = cell(1,2); LoR = cell(1,2); HiR = cell(1,2);
if ischar(varargin{1})
    [LD,HD,LR,HR] = wfilters(varargin{1}); 
    for k = 1:2
        LoD{k} = LD; HiD{k} = HD; LoR{k} = LR; HiR{k} = HR;
    end

elseif isstruct(varargin{1})
    if isfield(varargin{1},'w1') && isfield(varargin{1},'w2')
        for k = 1:2
            [LoD{k},HiD{k},LoR{k},HiR{k}] = ...
                wfilters(varargin{1}.(['w' int2str(k)]));
        end
    elseif isfield(varargin{1},'LoD') && isfield(varargin{1},'HiD') && ...
           isfield(varargin{1},'LoR') && isfield(varargin{1},'HiR')
        for k = 1:2
            LoD{k} = varargin{1}.LoD{k}; HiD{k} = varargin{1}.HiD{k};
            LoR{k} = varargin{1}.LoR{k}; HiR{k} = varargin{1}.HiR{k};
        end
    else
        error(message('Wavelet:FunctionArgVal:Invalid_ArgVal'));
    end
        
elseif iscell(varargin{1})
    if ischar(varargin{1}{1})
        for k = 1:2
            [LoD{k},HiD{k},LoR{k},HiR{k}] = wfilters(varargin{1}{k});
        end
    else
        LoD(1:end) = varargin{1}(1); HiD(1:end) = varargin{1}(2);
        LoR(1:end) = varargin{1}(3); HiR(1:end) = varargin{1}(4);
    end
else
    
end
nextArg = 2;

dwtEXTM = 'sym';
while nbIn>=nextArg
    argName = varargin{nextArg};
    argVal  = varargin{nextArg+1};
    nextArg = nextArg + 2;
    switch argName
        case 'mode' , dwtEXTM = argVal;
    end
end

% Initialization.
if isempty(X) , varargout{1} = []; return; end
sX = size(X);
X = double(X);
sizes = zeros(level+1,length(sX));
sizes(level+1,:) = sX;

for k=1:level
    dec = decFUNC(X,LoD,HiD,dwtEXTM);
    X = dec{1,1,1};
    sizes(level+1-k,:) = size(X);
    dec = reshape(dec,4,1,1);
    if k>1
        cfs(1) = [];
        cfs = cat(1,dec,cfs);
    else
        cfs = dec;
    end
end

WT.sizeINI = sX;
WT.level = level;
WT.filters.LoD = LoD;
WT.filters.HiD = HiD;
WT.filters.LoR = LoR;
WT.filters.HiR = HiR;
WT.mode = dwtEXTM;
WT.dec = cfs;
WT.sizes = sizes;
varargout{1} = WT;

%-------------------------------------------------------------------------%
function dec = decFUNC(X,LoD,HiD,dwtEXTM)

dec = cell(2,2);
permVect = [];
[a_Lo,d_Hi] = wdec1D(X,LoD{1},HiD{1},permVect,dwtEXTM);
permVect = [2,1,3];
[dec{1,1},dec{1,2}] = wdec1D(a_Lo,LoD{2},HiD{2},permVect,dwtEXTM);
[dec{2,1},dec{2,2}] = wdec1D(d_Hi,LoD{2},HiD{2},permVect,dwtEXTM);
%-------------------------------------------------------------------------%
function [L,H] = wdec1D(X,Lo,Hi,perm,dwtEXTM)

if ~isempty(perm) , X = permute(X,perm); end
sX = size(X);
if length(sX)<3 , sX(3) = 1; end
lf = length(Lo);
lx = sX(2);
lc = lx+lf-1;
switch dwtEXTM
    case 'zpd'             % Zero extension.
        
    case {'sym','symh'}    % Symmetric extension (half-point).
        X = [X(:,lf-1:-1:1,:) , X , X(:,end:-1:end-lf+1,:)];
        
    case 'sp0'             % Smooth extension of order 0.
        X = [X(:,ones(1,lf-1),:) , X , X(:,lx*ones(1,lf-1),:)];
        
    case {'sp1','spd'}     % Smooth extension of order 1.
        Z = zeros(sX(1),sX(2)+ 2*lf-2,sX(3));
        Z(:,lf:lf+lx-1,:) = X;
        last = sX(2)+lf-1;
        for k = 1:lf-1
            Z(:,last+k,:) = 2*Z(:,last+k-1,:)- Z(:,last+k-2,:);
            Z(:,lf-k,:)   = 2*Z(:,lf-k+1,:)- Z(:,lf-k+2,:);
        end
        X = Z; clear Z;
        
    case 'symw'            % Symmetric extension (whole-point).
        X = [X(:,lf:-1:2,:) , X , X(:,end-1:-1:end-lf,:)];
        
    case {'asym','asymh'}  % Antisymmetric extension (half-point).
        X = [-X(:,lf-1:-1:1,:) , X , -X(:,end:-1:end-lf+1,:)];        
        
    case 'asymw'           % Antisymmetric extension (whole-point).
        X = [-X(:,lf:-1:2,:) , X , -X(:,end-1:-1:end-lf,:)];

    case 'rndu'            % Uniformly randomized extension.
        X = [randn(sX(1),lf-1,sX(3)) , X , randn(sX(1),lf-1,sX(3))];        
                        
    case 'rndn'            % Normally randomized extension.
        X = [randn(sX(1),lf-1,sX(3)) , X , randn(sX(1),lf-1,sX(3))];        
                
    case 'ppd'             % Periodized extension (1).
        X = [X(:,end-lf+2:end,:) , X , X(:,1:lf-1,:)];
        
    case 'per'             % Periodized extension (2).
        if rem(lx,2) , X = [X , X(:,end,:)]; end
        X = [X(:,end-lf+2:end,:) , X , X(:,1:lf-1,:)];        
end
L = convn(X,Lo);
H = convn(X,Hi);
clear X
switch dwtEXTM
    case 'zpd'
    otherwise
        lenL = size(L,2);
        first = lf; last = lenL-lf+1;
        L = L(:,first:last,:); H = H(:,first:last,:);
        lenL = size(L,2);
        first = 1+floor((lenL-lc)/2);  last = first+lc-1;
        L = L(:,first:last,:); H = H(:,first:last,:);
end
if isequal(dwtEXTM,'per')
    first = 1; last = lx;
    L = L(:,first:last,:);
    H = H(:,first:last,:);
end

if ~isempty(perm)
    L = permute(L,perm);
    H = permute(H,perm);
end
%-------------------------------------------------------------------------%


