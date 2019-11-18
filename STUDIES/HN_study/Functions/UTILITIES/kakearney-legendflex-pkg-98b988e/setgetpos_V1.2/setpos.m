function setpos(h,fmt,href)
% SETPOS Set graphics object position in a flexible way.
%   SETPOS(H,FMT) sets the position property of graphics object 
%   with handle H, according to FMT that can be expressed using different
%   units. H must have a "Position' property.
%
%   FMT is a char array containing 4 strings separated by colon or space.
%   The format of each string is one of "%1c%f%2c" or "%1c%d%2c" where the 
%   first optional argument is "+" or "-", the second one is a number and 
%   the last one is two characters that specify the unit as :
%
%           px          for Pixels
%           nz          for Normalized
%           in          for Inches
%           cm          for Centimeters
%           pt          for Points
%           ch          for Characters
%           [] (empty)  for Current units [See get(H,'units')]
%
%   For better rendering, SETPOS can be included into the "Createfcn" or 
%   "Resizefcn" properties of the graphical object H.
%
%   Any string value of FMT can be replaced by a single '#' to keep the current 
%   value of the corresponding parameter.
%
%   The first optional argument of FMT is used to increase ('+') or 
%   decrease ('-') the corresponding value.
%
%   Note that SETPOS(H,FMT) works as set(H,'Position',FMT) when FMT is 
%   a 4 double values vector.
%
%   SETPOS(H,FMT,HREF) sets the position of the graphics object H according to 
%   FMT, but using the position of the graphics object HREF as reference instead 
%   of the parent of H. HREF must be a valid handle and must have a "Position" 
%   property (except for the Root object). Note that this should only affect 
%   Left&Bottom (1st&2nd) element of the "Position" vector of H.
%
%   See also GETPOS, SET, GET.

%   Author: Jérôme Briot, Matlab 6.1.0.450 (R12.1)
%   Contact: dutmatlab@yahoo.fr
%   Revision: 1.0 (12-Feb-2007)
%             1.1 (14-Feb-2007) Third input argument HREF added.
%                               Minor corrections in the help section.
%             1.2 (21-Feb-2007) Bug fixed if HREF is the Root object
%                               Examples removed from the help section
%   Comments:
%

% Check the number of input arguments

narginchk(2,3);


% Check if H is a graphics object handle
if ~ishandle(h)
    error('First argument must be a graphic object handle in SETPOS(H,FMT)');
end

% If FMT is a 4x1 double vector then SETPOS works as SET(H,'Position',FMT)
if isnumeric(fmt) & numel(fmt(:))==4
    
    set(h,'position',fmt)
    return

% If FMT is not a double vector, check if it's a char string
elseif ~ischar(fmt)

    error('FMT argument must be a string or a 4 elements vector in SETPOS(H,FMT)');
    
end

if nargin==2 % SETPOS(H,FMT)
    
    %HREF = parent of H
    href=get(h,'parent');
    
elseif nargin==3 % SETPOS(H,FMT,HREF)
    
    if ~ishandle(href) % Check if HREF is a valid handle
        error('HREF must be a valid handle of a graphics object in SETPOS(H,FMT,HREF)')
    end
    
end

flag_href=0;
% Don't use HREF position if it is the parent of H        
if href~=get(h,'parent') 
    flag_href=1;
end

% Extract 4 char strings from FMT
M=strread(fmt,'%s','delimiter',' ,','emptyvalue',0);

% Store the current unit of the graphics object H
current_unit=get(h,'units');
% Store the current unit of the reference object HREF
current_ref_unit=get(href,'units');

% List available units
available_units={'inches' 'centimeters' 'normalized' 'points' 'pixels' 'characters'};

flag=zeros(1,4);

% Decode elements of FMT
for n=1:numel(M)    
    
    % If FMT(n) is not a "#"
    if ~strcmp(M{n},'#')
        
        % Check if FMT(n) is +%... or -%...
        if strncmp(M{n},'+',1)
            flag(n)=1;
            M{n}(1)=[]; % Remove '+' char     
        elseif strncmp(M{n},'-',1)
            flag(n)=-1;
            M{n}(1)=[]; % Remove '-' char  
        end
        
        % Separate value and unit from FMT(n)
        [val(n),temp_unit]=strread(M{n},'%f%s');
        
        % If the unit is not specified in FMT(n)
		if isempty(temp_unit)
            
            unit{n}=current_unit; % Set the units to the current one
            
        % Else check if the units paramter is valid
 		else idx=strcmpi(temp_unit,{'in' 'cm' 'nz' 'pt' 'px' 'ch'});
            
            if ~any(idx)
                error('Units must be one of "in", "cm", "nz", "pt", "px" or "ch"')
            end
            
            unit{n}=available_units{idx}; % Set the units to one from the list
                
		end
        
    end           
                
end

% Set position of H using decoded FMT 
for n=1:numel(M)
    
    % If FMT(n) is not a "#" => value to modify
    if ~strcmp(M{n},'#')
        
        % Modify the "Units" property of H 
        set(h,'units',unit{n});
        % Modify the "Units" property of HREF
        set(href,'units',unit{n});
        % Get the current "Position" vector of H
        position_in_unit=get(h,'position');
        % Get the current "Position" vector of HREF
        if (isnumeric(href) && ~href) || (isgraphics(href) && isequal(href, groot)) % HREF is the Root object (no 'Position' property)
            position_ref_unit=get(href,'screensize'); %%% Should be safe here !
        else position_ref_unit=get(href,'position');
        end
        if ~flag % No "+" or "-"
            
            if any(n==[1 2])
                % If HREF is specified and is not the parent of H, flag_href=1 else flag_href=0
                position_in_unit(n)=val(n)+position_ref_unit(n)*flag_href;
            else position_in_unit(n)=val(n);
            end
            
        elseif any(n==[3 4]) % "+" or "-" and FMT(n) is width or height
        
            position_in_unit(n)=position_in_unit(n)+val(n)*flag(n);
            
        else % "+" or "-" and FMT(n) is left or bottom
            
            position_in_unit(n)=position_in_unit(n)+val(n)*flag(n);
            position_in_unit(n+2)=position_in_unit(n+2)-val(n)*flag(n);
            
        end
        
        % Modify the "Position" property of H
        set(h,'position',position_in_unit)
        
    end
        
end

% Restore the unit of the graphics object H
set(h,'units',current_unit);
% Restore the unit of the reference object HREF
set(href,'units',current_ref_unit);