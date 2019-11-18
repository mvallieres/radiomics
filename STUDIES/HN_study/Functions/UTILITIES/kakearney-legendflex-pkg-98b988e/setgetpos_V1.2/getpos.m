function [pos,unit]=getpos(h,fmt,href,opt)
% GETPOS Get graphics object position in a flexible way.
%   GETPOS(H,FMT) gets the position property of graphics object 
%   with handle H, according to FMT that can be expressed using different
%   units. H must have a "Position" property.
%
%   FMT is a char array containing four "%2c" strings separated by colon or
%   space. The two characters specify the unit as :
%
%           px  for Pixels
%           nz  for Normalized
%           in  for Inches
%           cm  for Centimeters
%           pt  for Points
%           ch  for Characters
%
%   If FMT is only one format string from the above list, all returned values are
%   expressed using this unit.
%
%   Any string value of FMT can be replaced by a single '#' to not retrieve the
%   corresponding value. The returned value is NaN except if the optional last 
%   argument OPT is set to "compact" in GETPOS(H,FMT,[HREF],OPT).
%
%   Note that GETPOS(H) works as get(H,'Position') and return the position 
%   vector in the current unit of the graphics object H.
%
%   GETPOS(H,FMT,HREF,['compact']) gets the position of the graphics object H according 
%   to FMT, but using the position of the graphics object HREF as reference instead 
%   of the parent of H. HREF must be a valid handle and must have a "Position" 
%   property (except for the Root object). Returned values may be negative or 0.
%
%   [POS,UNIT]=GETPOS(H,...) returns an additional output argument UNIT that 
%   contained the unit list of the output vector position POS. It may be safer 
%   when different units are used.
% 
%   See also SETPOS, SET, GET.

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

narginchk(1,4);


% Check if H is a graphics object handle
if ~ishandle(h)
    error('First argument must be a graphic object handle');
end

% Store the current unit of the graphics object H
current_unit=get(h,'units');

% Init variables
unit={current_unit current_unit current_unit current_unit};
pos=[nan nan nan nan];

% If FMT input argument is not specified, works as GET(H,'Position')
if nargin==1
    pos=get(h,'position');
    return
end

% Check if FMT is a char string
if ~ischar(fmt)
	error('Second argument must be a string in GETPOS(H,FMT)')
end  

if nargin==2 % GETPOS(H,FMT)
    
    href=get(h,'parent');
    opt='full';
    
elseif nargin==3
    
    if ishandle(href) % GETPOS(H,FMT,HREF)
        
        opt='full';
        
    elseif strcmpi(href,'compact') % GETPOS(H,FMT,"compact")
        
        href=get(h,'parent');
        opt='compact';
        
    else % GETPOS(H,FMT,???)
        error('Wrong third argument in GETPOS(H,FMT,???). Must be a valid handle or "compact"');
        
    end
    
elseif nargin==4 % GETPOS(H,FMT,HREF,OPT)
    
    if ~ishandle(href) 
        error('Third argument must be a valid handle in GETPOS(H,FMT,HREF,OPT)');
    end
        
    if ~strcmpi(opt,'compact') 
        error('Last argument must be "compact" in GETPOS(H,FMT,HREF,OPT)'); 
    end
    
end

flag_href=0;
% Don't use HREF position if it is the parent of H        
if href~=get(h,'parent')
    href=h;
    flag_href=1;
end

% Store the current unit of the reference object HREF
current_ref_unit=get(href,'units');

% Extract 4 char strings from FMT
M=strread(fmt,'%s','delimiter',' ,');

% Only one FMT requested for output
if numel(M)==1
    [M{2:4}]=deal(M{1});    
end

% List available units
available_units={'inches' 'centimeters' 'normalized' 'points' 'pixels' 'characters'};

% Decode elements of FMT
for n=1:numel(M) 
    
    % If FMT(n) is not a "#"
    if ~strcmp(M{n},'#')
        
        % Check if the units paramter is valid
        idx=strcmpi(M{n},{'in' 'cm' 'nz' 'pt' 'px' 'ch'});
        
        if ~any(idx)
            error('Units must be one of "in", "cm", "nz", "pt", "px" or "ch"')
        end
        
        unit{n}=available_units{idx}; % Set the units to one from the list 
            
    end
    
end

% Get position of H using decoded FMT 
for n=1:numel(M)    

    % If FMT(n) is not a "#" => get the value
    if ~strcmp(M{n},'#')
       
        % Modify the "Units" property of H 
        set(h,'units',unit{n});
        % Modify the "Units" property of HREF
        set(href,'units',unit{n});
        % Get the current "Position" vector of H
        temp=get(h,'position');
        % Get the current "Position" vector of HREF
        if strcmp(get(href, 'type'), 'root') % HREF is the Root object (no 'Position' property)
            temp_href=get(href,'screensize'); %%% Should be safe here !
        else temp_href=get(href,'position');
        end
        % Get and store the specified field from the "Position" vector
        % If HREF is specified and is not the parent of H, flag_href=1 else flag_href=0
        pos(n)=temp(n)-temp_href(n)*flag_href;
        
    end
    
end

% Check for compact output format 
if strcmpi(opt,'compact')
    pos(isnan(pos))=[];
end

% Restore the unit of the graphics object H
set(h,'units',current_unit);
% Restore the unit of the reference object HREF
set(href,'units',current_ref_unit);