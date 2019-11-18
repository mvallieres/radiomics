function X = indwt2(W,varargin)
%INDWT2 Inverse non-decimated 2-D wavelet transform.
%   INDWT2 performs a multilevel non-decimated 2-D wavelet reconstruction 
%   starting from a multilevel non-decimated 2-D wavelet decomposition.
%
%   In addition, you can use INDWT2 to simply extract coefficients 
%   from a multilevel non-decimated 2-D wavelet decomposition (see below).
%
%   C = INDWT2(W,TYPE,N) computes the reconstructed or the extracted
%   components at level N of a non-decimated 2-D wavelet decomposition. 
%   N must be a positive integer less or equal to the level of the
%   decomposition.
%   The valid values for TYPE are:
%       - A group of 2 chars 'xy', one per direction, with 'x' and 'y' 
%         in the set {'a','d','l','h'} or in the corresponding upper case  
%         set {'A','D','L','H'}), where 'A' (or 'L') stands for low pass 
%         filter and 'D' (or 'H') stands for high pass.
%       - The char 'd' (or 'h' or 'D' or 'H') gives directly the sum of 
%         the components different from the low pass one.
%   For extraction purpose, the valid values for TYPE are the same as
%   above prefixed by 'c' or 'C'.
%
%   C = INDWT2(W,TYPE) is equivalent to C = INDWT2(W,TYPE,N)
%   with N equal to the level of the decomposition.
%
%   X = INDWT2(W), X = INDWT2(W,'a',0) or X = INDWT(W,'ca',0)  
%   reconstructs the matrix X based on the non-decimated 2-D wavelet  
%   decomposition structure W.
%
%   Examples:
%       load noiswom
%       W = ndwt2(X,3,'db1');
%       A = cell(1,3);
%       for k=1:3, A{k} = indwt2(W,'aa',k); end
%       figure; colormap(pink(255))
%       subplot(2,3,2);image(X);
%       for k=1:3 
%           subplot(2,3,k+3);image(A{k});
%       end
%       D = indwt2(W,'d',1);
%       figure; colormap(pink(255));imagesc(abs(D))
%       A1 = indwt2(W,'aa',1); 
%       D1 = indwt2(W,'d',1);
%       E1 = X-A1-D1;
%       err1 = max(abs(E1(:)))
%       A2 = indwt2(W,'aa',2); 
%       D2 = indwt2(W,'d',2);
%       E2 = X-A2-D2;
%       err2 = max(abs(E2(:)))
%
%   See also DWTMODE, NDWT2, WAVEINFO.

%   M. Misiti, Y. Misiti, G. Oppenheim, J.M. Poggi 09-Dec-2008.
%   Last Revision: 20-Dec-2010.
%   Copyright 1995-2010 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $ $Date: 2011/04/12 20:43:22 $

% Check arguments.
nbIN = nargin-1;
idxCFS  = -1;
cfsFLAG = false;
if nbIN>0
    nbCELL = numel(W.dec);
    type = varargin{1};
    if ~ischar(type)
        error(message('Wavelet:FunctionArgVal:Invalid_ArgTyp'))
    end
    type = upper(type);
    cfsFLAG = isequal(upper(type(1)),'C');
    if cfsFLAG , type = type(2:end); end
    switch type
        case {'D','H'} ,           idxCFS = 0;
        case {'AA','LL','A','L'} , idxCFS = 1;
        case {'AD','LH'} ,         idxCFS = 2;
        case {'DA','HL'} ,         idxCFS = 3;
        case {'DD','HH'} ,         idxCFS = 4;
    end
    if nbIN>1 , levREC = varargin{2}; else levREC = W.level; end
        
    if idxCFS>1
        idxCFS = idxCFS + 3*(W.level-levREC);
        if ~cfsFLAG
            for j=1:nbCELL
                if ~isequal(j,idxCFS);
                    W.dec{j} = zeros(size(W.dec{j}));
                end
            end
        else
            X = W.dec{idxCFS};   % Coefficients
            return
        end
        
    elseif idxCFS==1   % Approximations (AA or LL)
        if cfsFLAG && levREC==W.level 
            X = W.dec{1}; 
            return; % Coefficients of Approximation at level MAX
        end
        idxMinToKill = 1 + 3*(W.level-levREC)+1;
        for j=idxMinToKill:nbCELL
            W.dec{j} = zeros(size(W.dec{j}));
        end
                
    elseif idxCFS==0
        idxMaxToKill = 1 + 3*(W.level-levREC);
        for j=1:idxMaxToKill
            W.dec{j} = zeros(size(W.dec{j}));
        end
        
    else
        
    end
end

% Initialization.
Lo  = W.filters.LoR;
Hi  = W.filters.HiR;
dwtEXTM = W.mode;
perFLAG = isequal(dwtEXTM,'per');
cfs   = W.dec;
sizes = W.sizes;
level = W.level;

maxloop = level;
if idxCFS==1 && cfsFLAG , maxloop = (level-levREC); end

idxBeg = 1;
for k=1:maxloop
    idxEnd = idxBeg+3;
    dec = reshape(cfs(idxBeg:idxEnd),2,2);
    sizerec = sizes(k+1,:);
    X   = recFUNC(dec,sizerec,Lo,Hi,perFLAG);
    cfs(1:idxEnd-1) = {[]};
    cfs{idxEnd} = X;
    idxBeg = idxEnd;
end

if abs(idxCFS)==1 && ~cfsFLAG && length(W.sizeINI)==3
    % X = uint8(X);
end
%-----------------------------------------------------------------------%
function X = recFUNC(dec,sINI,Lo,Hi,perFLAG)

% Reconstruction.
perm = [2,1,3];
W = cell(1,2);
for i = 1:2
    W{i} = wrec1D(dec{i,1},Lo{2},perm,perFLAG) + ...
        wrec1D(dec{i,2},Hi{2},perm,perFLAG);
end
X = (wrec1D(W{1},Lo{1},[],perFLAG) + wrec1D(W{2},Hi{1},[],perFLAG))/4;

% Extraction of central part
sREC = size(X);
F = floor((sREC-sINI)/2);
C = ceil((sREC-sINI)/2);
X = X(1+F(1):end-C(1),1+F(2):end-C(2),:);
%-----------------------------------------------------------------------%
function X = wrec1D(X,F,perm,perFLAG)

if ~isempty(perm) , X = permute(X,perm); end
if perFLAG
    nb = length(F)-1;
    X = [X X(:,1:nb,:)];
end
X = convn(X,F);
if ~isempty(perm) , X = permute(X,perm); end
%-----------------------------------------------------------------------%
