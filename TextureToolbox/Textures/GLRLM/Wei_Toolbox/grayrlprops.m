function stats = grayrlprops(varargin)

%GRAYCOPROPS Properties of gray-level run-length matrix.
%  -------------------------------------------
%  STATS = GRAYCOPROPS(GLRLM,PROPERTIES) Each element in  GLRLM, (r,c),
%   is the probability occurrence of pixel having gray level values r, run-length c in the image.
%   GRAYCOPROPS is to calculate PROPERTIES.
%  -------------------------------------------
%  Requirements:
%  -------------------------------------------
%   GLRLM mustbe an cell array of valid gray-level run-length
%   matrices.Recall that a valid glrlm must be logical or numerical.
%  -------------------------------------------
%  Current supported statistics include:
%  -------------------------------------------
%   Short Run Emphasis (SRE)
%   Long Run Emphasis (LRE)
%   Gray-Level Nonuniformity (GLN)
%   Run Length Nonuniformity (RLN)
%   Run Percentage (RP)
%   Low Gray-Level Run Emphasis (LGRE)
%   High Gray-Level Run Emphasis (HGRE)
%   Short Run Low Gray-Level Emphasis (SRLGE)
%   Short Run High Gray-Level Emphasis (SRHGE)
%   Long Run Low Gray-Level Emphasis (LRLGE)
%   Long Run High Gray-Level Emphasis (LRHGE)
%  --------------------------------------------
%  Reference:
%  --------------------------------------------
%   Xiaoou Tang,Texture Information in Run-Length Matrices
%   IEEE TRANSACTIONS ON IMAGE PROCESSING, VOL.7, NO.11,NOVEMBER 1998
% ---------------------------------------------
%  See also GRAYRLMATRIX.
% ---------------------------------------------
% Author:
% ---------------------------------------------
%    (C)Xunkai Wei <xunkai.wei@gmail.com>
%    Beijing Aeronautical Technology Research Center
%    Beijing %9203-12,10076
% ---------------------------------------------
% History:
% ---------------------------------------------
% Creation: beta         Date: 01/10/2007
% Revision: 1.0          Date: 12/11/2007
% 1.Accept cell input now
% 2.Using MATLAB file style
% 3.Fully vectorized programming
% 4.Fully support the IEEE reference
% 5. ...


% Check GLRLM
[GLRLM numGLRLM] = ParseInputs(varargin{:});

% Initialize output stats structure.
% 11 statistics for each GLRLM
numStats = 11;

% % count number of GLRLM
% numGLRLM = length(GLRLM);

% Initialization default 4*11 matrix
stats = zeros(numGLRLM,numStats);

for p = 1 : numGLRLM
    %N-D indexing not allowed for sparse.

    if numGLRLM ~= 1
        % transfer to double matrix
        tGLRLM = GLRLM{p};
    else
        tGLRLM = GLRLM;
    end
    %     if numGLRLM ~= 1
    %         % transfer to double matrix
    %         tGLRLM = normalizeGLRL(GLRLM{p});
    %     else
    %         tGLRLM = normalizeGLRL(GLRLM);
    %     end
    % Get row and column subscripts of GLRLM.  These subscripts correspond to the
    % pixel values in the GLRLM.
    s = size(tGLRLM);
    % colum indicator
    c_vector =1:s(1);
    % row indicator
    r_vector =1:s(2);
    % matrix element indicator
    % Matrix form col and row: using meshgrid, you should transpose before using
    % i.e. if tGLRLM is m*n, then this function return c_matrix n*m,
    % r_matrix n*m.
    [c_matrix,r_matrix] = meshgrid(c_vector,r_vector);

    % Total number of runs
    N_runs = sum(sum(tGLRLM));

    % total number of elements
    N_tGLRLM = s(1)*s(2);

    %--------------------Prepare four matrix for speedup--------------
    % 1.Gray Level Run-Length Pixel Number Matrix
    %     p_p = calculate_p_p(tGLRLM,c_matrix');

    % 2.Gray-Level Run-Number Vector
    %   This vector represents the sum distribution of the number of runs
    %   with gray level i.
    p_g = sum(tGLRLM);

    % 3.Run-Length Run-Number Vector
    %   This vector represents the sum distribution of the number of runs
    %   with run length j.
    p_r = sum(tGLRLM,2)';

    % 4.Gray-Level Run-Length-One Vector
    %
    % p_o = tGLRLM(:,1); % Not used yet
    % ----------------------End four matrix---------------------------
    %
    %------------------------Statistics-------------------------------
    % 1. Short Run Emphasis (SRE)
    SRE = sum(p_r./(c_vector.^2))/N_runs;
    % 2. Long Run Emphasis (LRE)
    LRE = sum(p_r.*(c_vector.^2))/N_runs;
    % 3. Gray-Level Nonuniformity (GLN)
    GLN = sum(p_g.^2)/N_runs;
    % 4. Run Length Nonuniformity (RLN)
    RLN = sum(p_r.^2)/N_runs;
    % 5. Run Percentage (RP)
    RP = N_runs/N_tGLRLM;
    % 6. Low Gray-Level Run Emphasis (LGRE)
    LGRE = sum(p_g./(r_vector.^2))/N_runs;
    % 7. High Gray-Level Run Emphasis (HGRE)
    HGRE = sum(p_g.*r_vector.^2)/N_runs;
    % 8. Short Run Low Gray-Level Emphasis (SRLGE)
    SGLGE =calculate_SGLGE(tGLRLM,r_matrix',c_matrix',N_runs);
    % 9. Short Run High Gray-Level Emphasis (SRHGE)
    SRHGE =calculate_SRHGE(tGLRLM,r_matrix',c_matrix',N_runs);
    % 10. Long Run Low Gray-Level Emphasis (LRLGE)
    LRLGE =calculate_LRLGE(tGLRLM,r_matrix',c_matrix',N_runs);
    % 11.Long Run High Gray-Level Emphasis (LRHGE
    LRHGE =calculate_LRHGE(tGLRLM,r_matrix',c_matrix',N_runs);
    %----------------insert statistics----------------------------
    stats(p,:)=[SRE LRE GLN RLN  RP LGRE HGRE SGLGE SRHGE LRLGE  LRHGE ];
end % end all run length matrixs

%   ----------------------Utility functions--------------------
%-----------------------------------------------------------------------------
% function glrl = normalizeGLRL(glrl)
%
% % Normalize glcm so that sum(glcm(:)) is one.
% if any(glrl(:))
%   glrl = glrl ./ sum(glrl(:));
% end
% function p_p = calculate_p_p(GLRLM,c) % Note: currently not used
%
% % p_p(i; j) = GLRLM(i,j)*j
% % Each element of the matrix represents the number of pixels of run length
% % j and gray-level i. Compared to the original matrix, the new matrix gives
% % equal emphasis to all length of runs in an image.
%
% term1 =  c; % j index in matrix size
% term2 = GLRLM;
% p_p = term1 .* term2;
%---------------------------------
function SGLGE =calculate_SGLGE(tGLRLM,r_matrix,c_matrix,N_runs)
% Short Run Low Gray-Level Emphasis (SRLGE):

term = tGLRLM./((r_matrix.*c_matrix).^2);
SGLGE= sum(sum(term))./N_runs;

%------------------------------------
function  SRHGE =calculate_SRHGE(tGLRLM,r_matrix,c_matrix,N_runs)
% Short Run High Gray-Level Emphasis (SRHGE):
%
term  = tGLRLM.*(r_matrix.^2)./(c_matrix.^2);
SRHGE = sum(sum(term))/N_runs;
%------------------------------------
function   LRLGE =calculate_LRLGE(tGLRLM,r_matrix,c_matrix,N_runs)
% Long Run Low Gray-Level Emphasis (LRLGE):
%
term  = tGLRLM.*(c_matrix.^2)./(r_matrix.^2);
LRLGE = sum(sum(term))/N_runs;
%---------------------------------------
function  LRHGE =calculate_LRHGE(tGLRLM,r_matrix,c_matrix,N_runs)
% Long Run High Gray-Level Emphasis (LRHGE):
%
term  = tGLRLM.*(c_matrix.^2).*(r_matrix.^2);
LRHGE = sum(sum(term))/N_runs;
%----------------------------------------

%-----------------------------------------------------------------------------
function [glrlm num_glrlm] = ParseInputs(varargin)
% check stability of inputs
%
% first receive all inputs
glrlm = varargin{:};
% get numbers total
num_glrlm=length(glrlm);
% then for each element, check its stability
for i=1:num_glrlm
    % The 'nonnan' and 'finite' attributes are not added to iptcheckinput because the
    % 'integer' attribute takes care of these requirements.
    % iptcheckinput(glrlm,{'cell'},{'real','nonnegative','integer'}, ...
    % mfilename,'GLRLM',1);
    iptcheckinput(glrlm{i},{'logical','numeric'},{'real','nonnegative','integer'},...
        mfilename,'GLRLM',1);
    % Cast GLRLM to double to avoid truncation by data type. Note that GLRLM is not an
    % image.
    if ~isa(glrlm,'double')
        glrlm{i}= double(glrlm{i});
    end
end
