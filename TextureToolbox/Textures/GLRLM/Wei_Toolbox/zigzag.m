function seq = zigzag(SI)
%
%  Description:
%  ------------
%  This function is used to build the corresponding sequences of a given
%  scaled gray level image matrix from 45' degree direction. The whole process is using zigzag method
%  It can handle nonsquare image matrix
%
% Author:
% -------
%    (C)Xunkai Wei <xunkai.wei@gmail.com>
%    Beijing Aeronautical Technology Research Center
%    Beijing %9203-12,10076
%
% History:
%  -------
% Creation: beta  Date: 01/11/2007
% Revision: 1.0   Date: 12/11/2007
% 
% Trick: all the sequence starts or ends lie on the boundary.

% initializing the variables
%----------------------------------
c = 1; % initialize colum indicator
r = 1; % initialize row   indicator

rmin = 1; % row   boundary checker
cmin = 1; % colum boundary checker

rmax = size(SI, 1); % get row numbers
cmax = size(SI, 2); % get colum numbers

%
i = 1; % counter for current ith element
j = 1; % indicator for determining sequence interval

% intialize sequence mark
sq_up_begin=1;

sq_down_begin=1;

% % Output contain value and its flag status
% the first row contain value
% the second row contain its flag
output = zeros(1, rmax * cmax);
% sequence counter
%
% % Position Matrix
% position =zeros(1, rmax * cmax);
%----------------------------------

while ((r <= rmax) && (c <= cmax))

    % for current point, judge its zigzag direction up 45, or down 45, or
    % 0,or down 90

    if (mod(c + r, 2) == 0)      % up 45 direction
        %  if we currently walk to the left first colum
        if (r == rmin)
            % First, record current point
            output(i) = SI(r, c);
            % if we walk to right last colum
            if (c == cmax)
                % add row number move straight down 90
                r   = r + 1;
                sq_up_end = i;
                sq_down_begin = i+1;
                seq{j}=output(sq_up_begin:sq_up_end);
                j = j + 1;
                %

            else
                % Continue to move to next (1,c+1) point
                % This next point should be the begin point of next sequence
                c = c + 1;
                sq_up_end = i;
                sq_down_begin = i+1;

                seq{j}=output(sq_up_begin:sq_up_end);

                j = j + 1;

            end;

            % add couter
            i = i + 1;
            % if we currently walk to the last column
        elseif ((c == cmax) && (r < rmax))
            % first record the point
            output(i) = SI(r, c);
            % then move straight down to next row
            r = r + 1;
            
            sq_up_end = i;
            seq{j}=output(sq_up_begin:sq_up_end);
            sq_down_begin =i+1;
            j=j+1;
                        
            % add counter
            i = i + 1;
            % all other cases i.e. nonboundary points
        elseif ((r > rmin) && (c < cmax))
            output(i) = SI(r, c);
            % move to next up 45 point
            r = r - 1;
            c = c + 1;
            % add counter
            i = i + 1;
        end;
        % down 45 direction
    else
        % if we walk to the last row
        if ((r == rmax) && (c <= cmax))
            % firstly record current point
            output(i) = SI(r, c);
            % move right to next point
            c = c + 1;
            sq_down_end = i;
            seq{j}=output(sq_down_begin:sq_down_end);
            sq_up_begin =i+1;
            j = j + 1;
            % add counter
            i = i + 1;
            % if we walk to the first column
        elseif (c == cmin)
            %first record current point
            output(i) = SI(r, c);
            %
            if (r == rmax)
                c = c + 1;
                
                sq_down_end = i;
                seq{j}=output(sq_down_begin:sq_down_end);
                sq_up_begin =i+1;
                j = j + 1;

            else
                r = r + 1;
                % record sequence end
                sq_down_end = i;
                seq{j}=output(sq_down_begin:sq_down_end);
                sq_up_begin =i+1;
                j = j + 1;

            end;

            i = i + 1;
            % all other cases without boundary point
        elseif ((r < rmax) && (c > cmin))
            %
            output(i) = SI(r, c);
            %             position(i) = sub2ind(SI,r,c);
            r = r + 1;
            c = c - 1;
            % keep down_info
            i = i + 1;
        end;

    end;

    if ((r == rmax) && (c == cmax))          % bottom right element
        output(i) = SI(r, c);
        sq_end = i;
        seq{j}=output(sq_end);
        %         position(i) = sub2ind(SI,r,c);
        break
    end;




end;