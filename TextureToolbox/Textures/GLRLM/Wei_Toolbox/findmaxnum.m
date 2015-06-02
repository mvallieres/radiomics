function maxnum=findmaxnum(seq)
%
%  this function is obtain the maximum numbers of the given sequence
%  note the sequence is stored in cell mode
%
%
% See also zigzag
%
% Author:
% ---------------------------------------------
%    (C)Xunkai Wei <xunkai.wei@gmail.com>
%    Beijing Aeronautical Technology Research Center
%    Beijing %9203-12,10076
% History:
% ---------------------------------------------
% Creation: beta  Date: 01/10/2007
% Revision: 1.0   Date: 12/11/2007
%
if iscell(seq)

    numseq=length(seq);
    maxnum=1;
    for i=1:numseq
        temp = seq{i};
        numseq = length(temp);
        if numseq > maxnum
            maxnum =numseq;
        end
    end
else
    error('I was only designed to handle cell sequence')
end