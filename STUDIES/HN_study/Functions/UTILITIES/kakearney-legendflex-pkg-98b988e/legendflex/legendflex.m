function varargout = legendflex(varargin)
%LEGENDFLEX Creates a more flexible legend
%
% legendflex(M, param1, val1, ...)
% legendflex(h, M, param1, val1, ...)
% [legend_h,object_h,plot_h,text_str] = legendflex(...)
% 
% This offers a more flexible version of the legend command.  It offers a
% different method of positioning the legend, as well as options to:
%
%   - organize legend text and symbols in a grid with a specified number of
%     rows and/or columns 
%   - rescale the horizontal space used by each legend symbol
%   - create multiple legends for the same axis
%   - add a title to the legend within the legend box
%
% Unlike in the default legend command, where the legend is positioned
% relative to the labeled objects' parent axis according to one of 16
% location strings, this function positions the legend based on two anchor
% points (one on either the figure or a child object of a figure, and one
% on the legend itself) and a buffer (or offset) between these two anchor
% points. The anchor points refer to the corners and centers of each
% side of the box surrounding the reference object and the legend itself;
% they can be refered to either as numbers (1-8, clockwise from northwest
% corner) or strings ('nw', 'n', 'ne', 'e', 'se', 's', 'sw', 'w').  The
% position of the legend is determined by these two points and the distance
% between them, defined in the 'buffer' variable, which by default is
% measured in pixels.  So the combination of
% 
%  (..., 'ref', gca, 'anchor', [3 3], 'buffer', [-10 -10])
%
% means that you want the northeast corner of the current axis to be
% aligned with the northeast corner of the legend, but with the legend
% shifted 10 pixels to the left and down. 
%
% This method of positioning can be particularly useful when labeling a
% figure that includes many subplots that share a common color scheme,
% where the "best" location for a legend is not necessarily within the
% bounds of an axis.  Unlike the legend command, the axes in the figure are
% never resized (and it is up to the user to check that the legend fits on
% the figure in the specified location).  In addition to being easier than
% manually positioning a legend, this function updates the legend location
% when the figure is resized, preserving the desired alignment.  The
% following anchor/buffer combinations, when used with the default
% reference and a buffer unit of pixels, approximately replicate the
% typical legend locations:
%
% Specifier              Anchor    Buffer
% 
% north                  [2 2]     [  0 -10]
% south                  [6 6]     [  0  10]
% east                   [4 4]     [-10   0]
% west                   [8 8]     [ 10   0]
% northeast              [3 3]     [-10 -10]
% northwest              [1 1]     [ 10 -10]
% southeast              [5 5]     [-10  10]
% southwest              [7 7]     [ 10  10]
% northoutside*          [2 6]     [  0  10]
% southoutside*          [6 2]     [  0 -10]
% eastoutside*           [3 8]     [ 10   0]
% westoutside*           [8 3]     [-10   0]
% northeastoutside*      [3 1]     [ 10   0]
% northwestoutside*      [1 3]     [-10   0]
% southeastoutside*      [5 7]     [ 10   0]
% southwestoutside*      [7 5]     [-10   0]  *placed outside axis rather
%                                              than resizing plot box
%
% This function should support all types of plot objects.
%
% Updates to labeled line and patch properties should be reflected in the
% legend.  In pre-R2014b versions of Matlab (those that use the old
% non-object graphics handles), properties of more complex legend labels,
% such as contours, quivers, bars, etc.) will also be synced to the legend;
% however, at this time, the code doesn't update properties for anything
% other than lines and patches in R2014b+ (haven't found a good way to
% listen for changes to the properties of the other graphics object types).
%
% A note on resizing: This function assigns a resize function to the parent
% figure to maintain the position of the legend (in terms of anchor
% location and buffer) as the figure size changes.  If you manually resize
% the legend, this function will respect changes to height, width, and
% units (though I don't recommend changing the units to 'normalized', as
% this can cause the text and symbols to overflow the legend box on
% resize).  It will not respect manual repositioning when resizing, since
% it assumes you want to maintain the anchor/buffer prescription used to
% create it.  Overall, I've tried to make this resize as unobtrusive as
% possible; if your figure already has a resize function at the time you
% apply it, that behavior is inherited, with the legend-resize called
% afterward.  If you plan to further modify the figure's resize function
% post-legendflex and want to maintain repositioning of the legends,
% retrieve the resize function via hfun = get(hfig, 'ResizeFcn'), pass it
% to the new resize function, and invoke it via feval(oldfun, h, ed), where
% h and ed are the default variables passed by a callback function.
%
% Input variables:
%
%   M:          cell array of strings, labels for legend
%
%   h:          handle of axis or handle(s) of object(s) to be labeled.  If
%               this is an axis handle, all children of the axis will be
%               included in the legend.  If not included, current axis is
%               used.
%
% Optional input variables (passed as parameter/value pairs): [default]
%
%   ncol:       number of columns, or 0 to indicate as many as necessary
%               given the # of labeled objects [1 if nrow is 0, 0
%               otherwise] 
%
%   nrow:       number of rows, or 0 to indicate as many as necessary
%               given the # of labeled objects [0]
%
%   ref:        handle of object used to position the legend. This can be
%               either a figure or a child object of a figure (and does not
%               need to relate in any way to the objects being labeled).
%               If not included, the reference will be to the axis that a
%               normal legend would be associated with (usually the parent
%               axis of the labeled objects, unless objects from multiple
%               axes are passed, in which case it's the parent object of
%               the first labeled object).
%
%   anchor:     1 x 2 array specifying which points of the reference object
%               and new legend, respectively, to anchor to each other.
%               Anchor points can be described using either numbers (in a 1
%               x 2 double array) or directional strings (in a 1 x 2 cell
%               array) as follows:
%               1:  'nw'    upper left corner
%               2:  'n'     center of top edge
%               3:  'ne'    upper right corner
%               4:  'e'     center of right edge
%               5:  'se'    bottom right corner
%               6:  's'     center of bottom edge
%               7:  'sw'    bottom left corner
%               8:  'w'     center of left edge         
%
%               [[3 3], i.e. {'ne' 'ne'}]
%
%   buffer:     1 x 2 array of horizontal and vertical distance,
%               respectively, from the reference anchor point to the legend
%               anchor point. Distance is measured in units specified by
%               bufferunit. [[-10 -10]]
%               
%   bufferunit: unit for buffer distance.  Note that this property only
%               affects the units used to position the legend, not the
%               units for the legend itself (which is always a fixed size,
%               based on the space needed to encapsulate the specified
%               symbols and text).  The 'normalized' units are normalized
%               to size of the figure. ['pixels']   
%
%   box:        'on' or 'off', specifies whether to enclose legend objects
%               in a box ['on']
%
%   xscale:     scalar value indicating scale factor to apply to the width
%               required by each symbol, relative to the size used by
%               legend. For example, 0.5 will shorten the lines/patches by
%               half. [1]
%
%   title:      A title string to be added inside the legend box, centered,
%               above all legend entries.  This can be either a string or a
%               cell array of strings; the latter will produce a multi-line
%               title. If empty, no title is added.  ['']
%
%   padding:    1 x 3 array, pixel spacing added to beginning of each
%               column (before symbol), between symbol and text, and after
%               text, respectively.  Usually, the default provides the
%               spacing typical of a regular legend, but occassionally the
%               extent properties wrap a little too close to text, making
%               things look crowded; in these cases you can try unsquishing
%               (or squishing, via use of negative values) things via this
%               parameter. [2 1 1]   
%
%   nolisten:   logical scalar.  If true, don't add the event listeners.
%               The event listeners update the legend objects when you
%               change a property of the labeled objects (such as line
%               style, color, etc.).  However, the updating requires the
%               legend to be redrawn, which can really slow things down,
%               especially if you're labelling lots of objects that get
%               changed together (if you change the line width of 100
%               labeled lines, the legend gets redrawn 100 times).  In more
%               recent releases, this also occurs when printing to file, so
%               I recommend setting this to true if you plan to print a
%               legend with a large number of labeled objects.  The legend
%               will still be redrawn on figure resize regardless of the
%               value of this parameter. [false]
%
%   In addition to these legendflex-specific parameters, this function will
%   accept any parameter accepted by the original legend function (e.g.
%   font properties) except 'location', 'boxon', 'boxoff', or 'hide'.
%
% Output variables:
%
%   legend_h:   handle of the legend axis.  It is not linked to an axis or
%               graphics objects in the same way as a Matlab legend.
%               However, on figure resize, all properties of the legend
%               objects are checked for changes, so adjusting the figure
%               size can re-link the legend to the labeled objects after
%               you have made changes to those objects.
% 
%   object_h:   handles of the line, patch, and text graphics objects
%               created in the legend 
% 
%   plot_h:     handles of the lines and other objects labeled in this
%               legend
% 
%   text_str:   cell array of the text strings used in the legend
%
%
% Example:
%
% % Replicating an example from legend.m:
%
% figure;
% b = bar(rand(10,5),'stacked'); colormap(summer); hold on
% x = plot(1:10,5*rand(10,1),'marker','square','markersize',12,...
%          'markeredgecolor','y','markerfacecolor',[.6 0 .6],...
%          'linestyle','-','color','r','linewidth',2); hold off
% lbl = {'Carrots','Peas','Peppers','Green Beans','Cucumbers','Eggplant'};
%
% % Rather than covering up data or resizing the axis, let's squeeze the
% % legend into the margin at the top of the figure;
%
% legendflex([b,x], lbl, 'ref', gcf, ...
%                        'anchor', {'n','n'}, ...
%                        'buffer',[0 0], ...
%                        'nrow',2, ...
%                        'fontsize',8);

% Copyright 2011-2014 Kelly Kearney

% Detemine whether HG2 is in use

hg2flag = ~verLessThan('matlab', '8.4.0');
r2016aflag = ~verLessThan('matlab', '9.0.0');

%-------------------
% Parse input
%-------------------
% 
% allinput = varargin; % Save for callback later
% 
% islegin = false(size(varargin));

% First inputs must be either:
% (M, ...)
% (h, M, ...)

narginchk(1,Inf);

% Split input into the variables that will be passed to legend (handles and
% labels) and everything else

handlepassed = all(ishandle(varargin{1})); % for HG1/HG2 

if handlepassed
    legin = varargin(1:2);
    if ~iscell(legin{2}) || ~all(cellfun(@ischar, legin{2}))
        error('Legend labels must be a cell array of strings');
    end
    pv = varargin(3:end);
else
    legin = varargin(1);
    if ~iscell(legin{1}) || ~all(cellfun(@ischar, legin{1}))
        if isnumeric(legin{1})
            error('Unable to parse input 1; check that handle(s) exist');
        else
            error('Legend labels must be a cell array of strings');
        end
    end
    pv = varargin(2:end);
end

% Parse my optional properties

if hg2flag
    defref = gobjects(0);
else
    defref = NaN;
end

p = inputParser;
p.addParameter('xscale',     1,         @(x) validateattributes(x, {'numeric'}, {'nonnegative','scalar'}));
p.addParameter('ncol',       0,         @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer'}));
p.addParameter('nrow',       0,         @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer'}));
p.addParameter('ref',        defref,    @(x) validateattributes(x, {'numeric','handle'}, {'scalar'}));
p.addParameter('anchor',     [3 3],     @(x) validateattributes(x, {'numeric','cell'}, {'size', [1 2]}));
p.addParameter('buffer',     [-10 -10], @(x) validateattributes(x, {'numeric'}, {'size', [1 2]}));
p.addParameter('bufferunit', 'pixels',  @(x) validateattributes(x, {'char'}, {}));
p.addParameter('box',        'on',      @(x) validateattributes(x, {'char'}, {}));
p.addParameter('title',      '',        @(x) validateattributes(x, {'char','cell'}, {}));
p.addParameter('padding',    [2 1 1],   @(x) validateattributes(x, {'numeric'}, {'size', [1 3]})); % 'nonnegative'
p.addParameter('nolisten',   false,     @(x) validateattributes(x, {'logical'}, {'scalar'}));

p.KeepUnmatched = true;

p.parse(pv{:});
Opt = p.Results;

% Any parameters that don't match mine are assumed to be a legend property.
%  If not, legend will handle the error when I call it.

Extra = p.Unmatched;
extra = [fieldnames(Extra) struct2cell(Extra)];
extra = extra';

% Validate that units and box inputs are correct

validatestring(Opt.bufferunit, {'pixels','normalized','inches','centimeters','points','characters'}, 'legendflex', 'bufferunit');
validatestring(Opt.box, {'on', 'off'}, 'legendflex', 'box');

% Translate anchor strings to numbers, if necessary

if iscell(Opt.anchor)
    [blah, Opt.anchor] = ismember(Opt.anchor, {'nw','n','ne','e','se','s','sw','w'});
    if ~all(blah)
        error('Anchor must be 1 x 2 cell array of strings: n, e, s, w, ne, nw, se, sw');
    end
else
    validateattributes(Opt.anchor, {'numeric'}, {'integer', '<=', 8}, 'legendflex', 'anchor');
end

% Create a temporary legend to get all the objects

S = warning('off', 'MATLAB:legend:PlotEmpty');
if r2016aflag
    % The new legend objects are pretty opaque... even diving into the 
    % undocumented properties, I haven't been able to find the handles of 
    % the legend sub-components (lines, text, etc).  So I need to stick to 
    % the legacy version, which creates an axis object rather than legend 
    % object. Legacy version has bug in text properties parsing, though, so 
    % need to work around that too: use the new-style legend object to get
    % proper text properties, then use those to alter the buggy old-style
    % legend.
    tmp = legend(legin{:}, extra{:}, 'location', 'northeast');
    textProps = {'FontAngle','FontName','FontSize','FontUnits','FontWeight','Interpreter'};
    tprop = get(tmp, textProps);
    delete(tmp);
    wtmp = warning('off', 'MATLAB:handle_graphics:exceptions:SceneNode'); % silence Latex interpreter thing
    [h.leg, h.obj, h.labeledobj, h.textstr] = legend(legin{:}, extra{:}, 'location', 'northeast');
    warning(wtmp);
    nobj = length(h.labeledobj);
    for it = 1:length(textProps)
        set(h.obj(1:nobj), textProps{it}, tprop{it});
    end
else
    [h.leg, h.obj, h.labeledobj, h.textstr] = legend(legin{:}, extra{:}, 'location', 'northeast');
    nobj = length(h.labeledobj);
end
warning(S);

if nobj == 0
    warning('Plot empty; no legend created');
    return
end

% There's a bug in R2014b-R2015a that causes rendering issues if a contour
% object is included in a legend and legend is called with more than one
% output. For some reason, the rendering issues disappear only if the
% contour object(s) is listed last in the legend.  So for now, my
% workaround for this is to change the order of the legend labels as
% necessary.  Issue appears to be fixed in 2015b.

iscont = strcmp(get(h.labeledobj, 'type'), 'contour');
cbugflag = ~verLessThan('matlab', '8.4.0') && verLessThan('matlab', '8.6.0') && any(iscont);

if cbugflag
    
    if length(legin) == 1
        legin = {h.labeledobj legin{1}};
    end
        
    delete(h.leg);
    
    [srt, isrt] = sort(iscont);
    legin{1} = legin{1}(isrt);
    legin{2} = legin{2}(isrt);
    
    [h.leg, h.obj, h.labeledobj, h.textstr] = legend(legin{:}, extra{:}, 'location', 'northeast');

end

% # rows and columns

if (Opt.ncol == 0) && (Opt.nrow == 0)
    Opt.ncol = 1;
    Opt.nrow = nobj;
elseif (Opt.ncol == 0)
    Opt.ncol = ceil(nobj./Opt.nrow);
elseif (Opt.nrow == 0)
    Opt.nrow = ceil(nobj./Opt.ncol);
end
if Opt.ncol*Opt.nrow < nobj
    error('Number of legend entries greater than specified grid allows; change ncol and/or nrow');
end

% Reference object

if hg2flag

    if isempty(Opt.ref)
    
        if all(ishandle(legin{1}))
            tmp = ancestor(legin{1}, 'axes');
            if iscell(tmp)
                Opt.ref = tmp{1}; 
            else
                Opt.ref = tmp(1);
            end
        else
            Opt.ref = gca;
        end

    end
else
    if isnan(Opt.ref)
        tmp = get(h.leg, 'UserData');
        Opt.ref = tmp.PlotHandle; 
    end
end
if ~ishandle(Opt.ref)
    error('Input ref must be a graphics handle');
end

% Box

Opt.box = strcmpi('on', Opt.box);

% Convert units to getpos abbreviations

unittable = {...
    'px'  'Pixels'
    'nz'  'Normalized'
    'in'  'Inches'
    'cm'  'Centimeters'
    'pt'  'Points'
    'ch'  'Characters'};
Opt.bufunit = unittable{strcmpi(unittable(:,2),Opt.bufferunit),1};

% Check for title

addtitle = ~isempty(Opt.title);

%-------------------
% New placement of
% everything in 
% legend
%-------------------

% Determine parent figure

figh = ancestor(Opt.ref, 'figure');
currax = get(figh, 'currentaxes'); 

% Calculate row height

legpospx = getpos(h.leg, 'px');

% rowHeight = legpospx(4)/nobj;
vmarginNm =  0.275/nobj;
vmarginPx = legpospx(4) * vmarginNm;

rowHeightNm = (1 - vmarginNm)/nobj;
rowHeight = rowHeightNm .* legpospx(4);

% Determine width needed for each text string

if nobj == 1
    textExtent = get(h.obj(1:nobj), 'Extent');
else
    textExtent = cell2mat(get(h.obj(1:nobj), 'Extent'));
end
textWidthPx  = textExtent(:,3) .* legpospx(3);
textHeightPx = textExtent(:,4) .* legpospx(4);
textWidthNm = textExtent(:,3);

% Calculate horizontal space needed for symbols

symbolWidthPx = textExtent(1,1) .* legpospx(3) * Opt.xscale;
symbolWidthNm = textExtent(1,1);

% Calculate column width needed for 2px-symbol-1px-text-1px

colWidth = zeros(Opt.ncol*Opt.nrow,1);
colWidth(1:nobj) = textWidthPx + symbolWidthPx + sum(Opt.padding);
colWidth = reshape(colWidth, Opt.nrow, Opt.ncol);
colWidth = max(colWidth,[],1);

% If title is added, figure out how much space it will need

if addtitle
    textProps = {'FontAngle','FontName','FontSize','FontUnits','FontWeight','Interpreter'};
    textVals = get(h.obj(1), textProps);
    ttlprops = [textProps; textVals];
    
    fpos = getpos(figh, 'px');
    figtmp = figure('units','pixels','position',[0 0 fpos(3:4)],'visible','off');
    axes('parent',figtmp,'position',[0 0 1 1],'xlim',[0 fpos(3)],'ylim',[0 fpos(4)]);
    tmp = text(0,0,Opt.title, ttlprops{:}, 'horiz', 'left', 'vert', 'bottom');
    ttlex = get(tmp, 'extent');
    ttlwidth = ceil(ttlex(3)) + 4; % Add a little padding
    ttlheight = ceil(ttlex(4));
    
    if ttlwidth > sum(colWidth)
        colWidth(end) = colWidth(end) + (ttlwidth-sum(colWidth));
    end
    close(figtmp);
end

% Locate bottom left corner of each legend symbol, text box, and title

xsymbnew = [0 cumsum(colWidth(1:end-1))]+Opt.padding(1);
ysymbnew = (rowHeight*Opt.nrow + vmarginPx)-(1:Opt.nrow)*rowHeight;
[xsymbnew, ysymbnew] = meshgrid(xsymbnew, ysymbnew);
xsymbnew = xsymbnew(1:nobj);
ysymbnew = ysymbnew(1:nobj);

xtext = xsymbnew + Opt.padding(2) + symbolWidthPx;
ytext = ysymbnew;% + 1;

xsymbold = zeros(nobj,1);
ysymbold = 1 - (1/nobj)*(1:nobj);

wnewleg = sum(colWidth);
hnewleg = rowHeight*Opt.nrow + vmarginPx;

if addtitle
    xttl = wnewleg/2;
    yttl = hnewleg;
    hnewleg = hnewleg + ttlheight;
end
    
% Get legend position in bufferunit and translate to pixels

legpos = positionleg(Opt.ref, wnewleg, hnewleg, Opt.anchor, Opt.buffer, Opt.bufunit);
tmpax = axes('units', Opt.bufferunit, 'position', legpos,'visible','off');
legpos = getpos(tmpax, 'px');
delete(tmpax);

%-------------------
% Create legend
%-------------------

% Create the legend axis

hnew.leg = axes('units', 'pixels', ...
               'position', legpos, ...
               'xlim', [0 legpos(3)], ...
               'ylim', [0 legpos(4)], ...
               'xtick', [], ...
               'ytick', [], ...
               'box', 'on', ...
               'parent', figh);

% Copy the text strings to the new legend
           
textProps = {'FontAngle','FontName','FontSize','FontUnits','FontWeight','Interpreter','HorizontalAlignment','VerticalAlignment'};
textVals = get(h.obj(1:nobj), textProps);

if hg2flag
    hnew.obj = gobjects(size(h.obj));
else
    hnew.obj = zeros(size(h.obj));
end
for it = 1:nobj
    props = [textProps; textVals(it,:)];
    hnew.obj(it) = text(xtext(it), ytext(it), h.textstr{it}, props{:}, ...
                        'horizontalalignment', 'left', ...
                        'verticalalignment', 'bottom');
end

% Copy the symbols to the new legend

nsymbol = length(h.obj) - nobj;

for ii = 1:nsymbol
    
    if strcmp(get(h.obj(nobj+ii), 'type'), 'hggroup')
        
        tag = get(h.obj(nobj+ii),'Tag');
        if ~isempty(tag)
            [blah, idx] = ismember(tag,h.textstr);
        end
        
        chld = findall(h.obj(nobj+ii), 'type', 'line', '-or', 'type', 'patch');       
        for ic = 1:length(chld)
            xy = get(chld(ic), {'xdata', 'ydata'});
            
            xnorm = xy{1}./symbolWidthNm;
            ynorm = (xy{2}- (1-idx*rowHeightNm))./rowHeightNm;

            xnew = xnorm * symbolWidthPx + xsymbnew(idx);
            ynew = ynorm * rowHeight     + ysymbnew(idx);
            
            set(chld(ic), 'xdata', xnew, 'ydata', ynew);
        end
        
        hnew.obj(nobj+ii) = copyobj(h.obj(nobj+ii), hnew.leg);
        
    else   
        
        hnew.obj(nobj+ii) = copyobj(h.obj(nobj+ii), hnew.leg);
        
        tag = get(h.obj(nobj+ii),'Tag');
        if ~isempty(tag) % assumes empty tags indicate repetition of previous tag (true pre-2014b)
            [blah, idx] = ismember(tag,h.textstr);
        end
        
        xy = get(h.obj(nobj+ii), {'xdata', 'ydata'});

        xnorm = xy{1}./symbolWidthNm;
        ynorm = (xy{2}- (1-idx*rowHeightNm))./rowHeightNm;

        xnew = xnorm * symbolWidthPx + xsymbnew(idx);
        ynew = ynorm * rowHeight     + ysymbnew(idx);

        set(hnew.obj(nobj+ii), 'xdata', xnew, 'ydata', ynew);
 
    end
    
end

% Add title

if addtitle
    text(xttl, yttl, Opt.title, ttlprops{:}, 'horiz', 'center', 'vert', 'bottom');
end

% Add box or hide axis

if Opt.box
    set(hnew.leg, 'box', 'on');
else
    set(hnew.leg, 'visible', 'off');
end

% Delete the temporary legend

delete(h.leg);

% Return focus to previously-current axis

set(figh, 'currentaxes', currax);
drawnow; % Not sure why this is necessary for the currentaxes to take effect, but it is

% Fix for vertical-alignment issue: This solution still isn't perfect, but
% it seems to help for most Interpreter-none and Interpreter-latex cases.
% The TeX interpreter still places sub- and superscripts too high/low... no
% robust fix found for that yet.
%
% TODO: need to add proper calcs for when title included
%
% Thanks to Søren Enemark for this suggestion.

if ~addtitle
    try % TODO: Crashing on some edge cases
        textobj = hnew.obj(1:nobj);
        yheight = get(hnew.leg, 'ylim');
        yheight = yheight(2);

        ylo = get(textobj(Opt.nrow), 'extent');
        ylo = ylo(2);
        yhi = get(textobj(1), 'extent');
        yhi = sum(yhi([2 4]));
        dy = yheight/2 - 0.5*(ylo + yhi);
        for ii = 1:length(textobj)
            pos = get(textobj(ii), 'position');
            set(textobj(ii), 'position', pos + [0 dy 0]);
        end
    end
end

%-------------------
% Callbacks and 
% listeners
%-------------------

% Save some relevant variables in the new legend axis's application data

Lf.ref        = Opt.ref;
Lf.w          = wnewleg;
Lf.h          = hnewleg;
Lf.anchor     = Opt.anchor;
Lf.buffer     = Opt.buffer;
Lf.bufunit    = Opt.bufunit;
Lf.bufferunit = Opt.bufferunit;
Lf.plotobj    = h.labeledobj;
Lf.legobj     = hnew.obj;

setappdata(hnew.leg, 'legflex', Lf);

% Resize listeners

addlistener(hnew.leg, 'Position', 'PostSet', @(src,evt) updatelegappdata(src,evt,hnew.leg));
if hg2flag && strcmp(Lf.ref.Type, 'figure')
    addlistener(Lf.ref, 'SizeChanged', @(src,evt) updatelegpos(src,evt,hnew.leg));
else
    addlistener(Lf.ref, 'Position', 'PostSet', @(src,evt) updatelegpos(src,evt,hnew.leg));
end
rsz = get(figh, 'ResizeFcn'); 
if isempty(rsz) % No previous resize function
    set(figh, 'ResizeFcn', @updatelegfigresize);
else 
    if ~iscell(rsz)
        rsz = {rsz};
    end
    hasprev = cellfun(@(x) isequal(x, @updatelegfigresize), rsz);
    if ~hasprev
        rsz = {rsz{:} @updatelegfigresize};
        set(figh, 'ResizeFcn', {@wrapper, rsz});
    end
end

if ~Opt.nolisten

    % Run the resync function if anything changes with the labeled objects

    objwatch = findall(h.labeledobj, 'type', 'line', '-or', 'type', 'patch');

    for ii = 1:length(objwatch)
        switch lower(get(objwatch(ii), 'type'))
            case 'line'
                triggerprops = {'Color','LineStyle','LineWidth','Marker','MarkerSize','MarkerEdgeColor','MarkerFaceColor'};
                addlistener(objwatch(ii), triggerprops, 'PostSet', @(h,ed) resyncprops(h,ed,hnew.leg));
            case 'patch'
                triggerprops = {'CData','CDataMapping','EdgeAlpha','EdgeColor','FaceAlpha','FaceColor','LineStyle','LineWidth','Marker','MarkerEdgeColor','MarkerFaceColor','MarkerSize'};
                addlistener(objwatch(ii), triggerprops, 'PostSet', @(h,ed) resyncprops(h,ed,hnew.leg));
        end
    end

end

    
%-------------------
% Output
%-------------------

out = {hnew.leg, hnew.obj, h.labeledobj, h.textstr};
varargout = out(1:nargout);


%***** Subfunctions *****

%------------------------
% Position new legend
%------------------------

function legpos = positionleg(href, w, h, anchor, buffer, bufunit)
% ap: position vector for reference object
% lp: position vector for legend

if strcmp(get(href, 'type'), 'figure')
    tmp = axes('parent', href,'position', [0 0 1 1],'visible','off');
    pos = getpos(tmp, bufunit);
    delete(tmp);
else
    pos = getpos(href, bufunit);
end

htmp = axes('units', 'pixels', 'position', [0 0 w h], 'visible','off');
lpos = getpos(htmp, bufunit);
delete(htmp);
w = lpos(3);
h = lpos(4);

% Find anchor locations on reference object

refxy = [...
    pos(1)          pos(2)+pos(4)
    pos(1)+pos(3)/2 pos(2)+pos(4)
    pos(1)+pos(3)   pos(2)+pos(4)
    pos(1)+pos(3)   pos(2)+pos(4)/2
    pos(1)+pos(3)   pos(2)
    pos(1)+pos(3)/2 pos(2)
    pos(1)          pos(2)
    pos(1)          pos(2)+pos(4)/2];

% How bottom left relates to each anchor point

shift = [...
    0       -h
    -w/2    -h
    -w      -h
    -w      -h/2
    -w      0
    -w/2    0
    0       0
    0       -h/2];

% Legend location

corner = refxy(anchor(1),:) + buffer + shift(anchor(2),:);
legpos = [corner w h];

%------------------------
% Listener functions
%------------------------

% If user manually resizes the legend, update the app data

function updatelegappdata(src, evt, legax)
if ishandle(legax)
    Lf = getappdata(legax, 'legflex');
    pos = getpos(legax, 'px');
    Lf.w = pos(3);
    Lf.h = pos(4);
    setappdata(legax, 'legflex', Lf);
end
% If reference object moves or resizes, reposition the legend appropriately

function updatelegpos(src, evt, legax)
if ishandle(legax) 
    Lf = getappdata(legax, 'legflex');
    legpos = positionleg(Lf.ref, Lf.w, Lf.h, Lf.anchor, Lf.buffer, Lf.bufunit);
    set(legax, 'Units', Lf.bufferunit, 'Position', legpos);
end

% Since figure resize can change axis size without actually triggering a
% listener, force this

function updatelegfigresize(src, evt)

allax = findall(src, 'type', 'axes');
for ii = 1:length(allax)
    isleg = ~isempty(getappdata(allax(ii), 'legflex'));
    if ~isleg
        pos = get(allax(ii), 'Position');
        set(allax(ii), 'Position', pos); % No change, just trigger PostSet
    end
end

% If plotted object changes, resync with legend

function resyncprops(src, evt, legax)

if ishandle(legax) % In case it's been deleted
    
    Lf = getappdata(legax, 'legflex');

    str = cellstr(num2str((1:length(Lf.plotobj))'));
    [htmp.leg, htmp.obj, htmp.labeledobj, htmp.textstr] = legend(Lf.plotobj, str);

    objtype = get(Lf.legobj, 'type');
    isline = strcmp(objtype, 'line');
    ispatch = strcmp(objtype, 'patch');
    ishg = strcmp(objtype, 'hggroup');
    hgidx = find(ishg);

    lobj = [Lf.legobj(isline) htmp.obj(isline)];
    pobj = [Lf.legobj(ispatch) htmp.obj(ispatch)];

    if ~isempty(hgidx)
        for ih = hgidx
            chldln1 = findall(Lf.legobj(ih), 'type', 'line');
            chldln2 = findall(htmp.obj(ih), 'type', 'line'); 

            lobj = [lobj; [chldln1 chldln2]];

            chldpa1 = findall(Lf.legobj(ih), 'type', 'patch');
            chldpa2 = findall(htmp.obj(ih), 'type', 'patch'); 

            pobj = [pobj; [chldpa1 chldpa2]];

        end
    end

    lprops = {'color','linestyle','linewidth','marker','markersize','markeredgecolor','markerfacecolor'};
    for il = 1:size(lobj,1)
        lvals = get(lobj(il,2), lprops);
        pv = [lprops; lvals];
        set(lobj(il,1), pv{:});
    end

    pprops = {'cdata','cdatamapping','edgealpha','edgecolor','facealpha','facecolor','linestyle','linewidth','marker','markeredgecolor','markerfacecolor','markersize'};
    for ip = 1:size(pobj,1)
        pvals = get(pobj(ip,2), pprops);
        pv = [pprops; pvals];
        set(pobj(ip,1), pv{:});
    end

    cmap = colormap(htmp.leg);
    colormap(legax, cmap);

    delete(htmp.leg);
end

% Wrapper to add multiple callback functions to resize

function wrapper(ObjH, EventData, fcnList)
for ii = 1:length(fcnList)
    feval(fcnList{ii}, ObjH, EventData);
end






    

