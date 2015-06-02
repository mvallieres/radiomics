function oneglrlm = rle_45(seq,NL)
% RLE   image gray level Run Length matrix for 45 and 135
% This file is to handle the zigzag scanned sequence for 45 or 135 degree
% direction. Note for 135, just swap the left and the right colum
% Author:
% ---------------------------------------------
%    (C)Xunkai Wei <xunkai.wei@gmail.com>
%    Beijing Aeronautical Technology Research Center
%    Beijing %9203-12,10076
% History:
%  -------
% Creation: beta  Date: 01/11/2007 
% Revision: 1.0   Date: 10/11/2007


% Assure row number is exactly the gray leve;
% number of seqence
m =length(seq);
% number to store the possible max coloums
n = findmaxnum(seq);
%

oneglrlm=zeros(NL,n);

for i=1:m
    x=seq{i};
    % run length Encode of each vector
    index = [ find(x(1:end-1) ~= x(2:end)), length(x) ];
    len = diff([ 0 index ]); % run lengths
    val = x(index);          % run values
    temp =accumarray([val;len]',1,[NL n]);% compute current numbers (or contribution) for each bin in GLRLM
    oneglrlm = temp + oneglrlm; % accumulate each contribution
end