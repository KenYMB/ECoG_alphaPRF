function hh = errorshade(varargin)
%ERRORSHADE Plot error shade along curve
%   ERRORSHADE(Y,E) plots Y and draws a vertical error shade at each element of
%   Y. The error shade is a distance of E(i) above and below the curve so
%   that the length of each side of shade is symmetric and 2*E(i) long.
%
%   ERRORSHADE(X,Y,E) plots Y versus X with symmetric vertical error shade
%   2*E(i) long. X, Y, E must be the same size. When they are vectors, each
%   error shade is a distance of E(i) above and below the point defined by
%   (X(i),Y(i)). When they are matrices, each error shade is a distance of
%   E(i,j) above and below the point defined by (X(i,j),Y(i,j)).
%
%   ERRORSHADE(X,Y,NEG,POS) plots X versus Y with vertical error shade
%   NEG(i)+POS(i) long specifying the lower and upper error shade. X and Y
%   must be the same size. NEG and POS must be the same size as Y or empty.
%   When they are vectors, each error shade is a distance of NEG(i) below and
%   POS(i) above the point defined by (X(i),Y(i)). When they are matrices,
%   each error shade is a distance of NEG(i,j) below and POS(i,j) above the
%   point defined by (X(i,j),Y(i,j)). When they are empty the error shade is
%   not drawn.
%
%   ERRORSHADE( ___ ,C) determines the colors for plot. The color is
%   applied to both plots and error shades. If C is a string, each plot is
%   filled with 'color'. 'color' can be 'r','g','b','c','m','y', 'w', 'k',
%   or a hexadecimal color code. If C is an array of strings it specifies
%   the color of the plots by indexing into the colormap. A 1-by-3 vector
%   or N-by-3 matrix C is always assumed to be RGB triplets specifying
%   color directlies.
%
%   ERRORSHADE( ___ ,Orientation) specifies the orientation of the error
%   shade. Orientation can be 'horizontal' or 'vertical'. When the
%   orientation is omitted the default is 'vertical'.
%
%   ERRORSHADE( ___ ,LineSpec) specifies the color, line style, and marker.
%   The color is applied to the data line and error shade. The line style
%   and marker are applied to the data line only.
%
%   ERRORSHADE(AX, ___ ) plots into the axes specified by AX instead of the
%   current axes.
%
%   H = ERRORSHADE( ___ ) returns handles to the Line and Patch objects
%   created. ERRORSHADE creates one object for vector input arguments and one
%   object per column for matrix input arguments.
%
%   Example: Draws symmetric error shade of unit standard deviation.
%      x = 1:10;
%      y = sin(x);
%      e = std(y)*ones(size(x));
%      errorshade(x,y,e)
% 
%   See also ERRORBAR.

%   20231030 Yuasa: coded based on ERRORBAR


% Look for a parent among the input arguments
[~, cax, args] = parseplotapi(varargin{:},'-mfilename',mfilename);

if isempty(cax) && ~isempty(args) && isa(args{1},'matlab.graphics.Graphics')
    % axescheck will not parse out first arg parents that are not axes,
    % e.g. charts. Manually extract these here so that they are not
    % confused for x/y data. Note that cax will be set back to empty below
    % by cax = ancestor(cax,'axes'); and en error will be thrown if the
    % parent is invalid.
    cax=args{1};
    args(1)=[];
end

% Errorshade requires at least two inputs
narginchk(2, inf);

% Separate Name/Value pairs from data inputs, convert LineSpec to
% Name/Value pairs, and filter out the orientation flag.
[pvpairs,args,nargs,msg,orientation,color] = parseargs(args);
if ~isempty(msg), error(msg); end
pvpairs = matlab.graphics.internal.convertStringToCharArgs(pvpairs);
% Check that we have the correct number of data input arguments.
if nargs < 2
    error(message('MATLAB:narginchk:notEnoughInputs'));
elseif nargs > 5 || ( nargs > 4 && ~isempty(args{end}) )
    error(message('MATLAB:narginchk:tooManyInputs'));
end

% Make sure all the data input arguments are real numeric data or a
% recognized non-numeric data type.
args = matlab.graphics.chart.internal.getRealData(args, true);

% Grab the X data if present.
if nargs >= 3
    % errorshade(x,y,e,...)
    x = args{1};
    args = args(2:end);
else
    % errorshade(y,e)
    x = [];
end

% Grab the Y data
y = checkSingleInput(args{1}, [], 'y');
sz = size(y);
n = sz(2);

% Now that we have the size of the YData, validate the size of the XData
x = checkSingleInput(x, sz, 'x');

% Grab the first delta inputs.
neg = args{2};

% Grab the second delta inputs.
if numel(args) >= 3
    % errorshade(x,y,neg,pos,...)
    pos = args{3};
    negArgumentName = 'neg';
    posArgumentName = 'pos';
else
    % errorshade(y,e) or
    % errorshade(x,y,e)
    pos = neg;
    negArgumentName = 'err';
    posArgumentName = 'err';
end

% Grab the remaining delta inputs and validate all data inputs.
switch orientation
    case 'both'
        if ~issorted(x,1,'monotonic')
            if issorted(y,1,'monotonic')
                orientation = 'horizontal';
            elseif sum(issortedarray(x,1,'monotonic')) < sum(issortedarray(y,1,'monotonic'))
                orientation = 'horizontal';
            elseif sum(issortedarray(x,1,'monotonic')) == sum(issortedarray(y,1,'monotonic'))
                sortedratex = sum(diff(x,1,1)>0,'all','omitnan') ./ sum(diff(x,1,1)<0,'all','omitnan');
                if sortedratex<=0, sortedratex = 1./sortedratex; end %#ok<BDSCI> 
                sortedratey = sum(diff(y,1,1)>0,'all','omitnan') ./ sum(diff(y,1,1)<0,'all','omitnan');
                if sortedratey<=0, sortedratey = 1./sortedratey; end %#ok<BDSCI> 
                if sortedratex < sortedratey
                    orientation = 'horizontal';
                end
            end
        end
end
switch orientation
    % errorshade(y,e,orientation) or
    % errorshade(x,y,e,orientation) or
    % errorshade(x,y,neg,pos,orientation)
    case 'horizontal'
        xneg = checkSingleInput(neg, sz, negArgumentName, x);
        xpos = checkSingleInput(pos, sz, posArgumentName, x);
        yneg = defaultEmptyDelta(y);
        ypos = defaultEmptyDelta(y);
    otherwise
        % Default to vertical if orientation isn't specified.
        xneg = defaultEmptyDelta(x);
        xpos = defaultEmptyDelta(x);
        yneg = checkSingleInput(neg, sz, negArgumentName, y);
        ypos = checkSingleInput(pos, sz, posArgumentName, y);
end

% Handle vectorized data sources and display names (back up obsolete variable)
extrapairs = cell(n,0);

% Prepare the parent for plotting.
if isempty(cax) || ishghandle(cax,'axes')
    showInteractionInfoPanel = isempty(cax) && isempty(get(groot,'CurrentFigure'));
    cax = newplot(cax);
    parax = cax;
    if showInteractionInfoPanel
        % Maybe open the Interaction Info Panel
        matlab.graphics.internal.InteractionInfoPanel.maybeShow(cax);
    end    
    hold_state = any(strcmp(cax.NextPlot,{'replacechildren','add'}));
else
    parax = cax;
    cax = ancestor(cax,'axes');
    hold_state = true;
end

% Configure the axes for non-numeric data.
xIsNumeric = isnumeric(x);
yIsNumeric = isnumeric(y);
matlab.graphics.internal.configureAxes(cax,x,y);

% Determine the Color and LineStyle property names
% If the Color/LineStyle is not specified use the _I property names so that
% the ColorMode or LineStyleMode properties are not toggled.
colorPropName = 'Color';
multiColor = ~isempty(color);
autoColor = ~multiColor & ~any(strcmpi('color',pvpairs(1:2:end)));
if autoColor
    colorPropName = 'Color_I';
end
stylePropName = 'LineStyle';
autoStyle = ~any(strcmpi('linestyle',pvpairs(1:2:end)));
if autoStyle
    stylePropName = 'LineStyle_I';
end

h = gobjects(2,n);
% Create the Line objects
for k = 1:n
    % extract data from vectorizing over columns
    xdata = {'XData', getColumn(x,k,xIsNumeric)};
    ydata = {'YData', getColumn(y,k,yIsNumeric)};
    
    stylepv={};
    if ~isempty(cax)
        [ls,c,m] = matlab.graphics.chart.internal.nextstyle(cax,autoColor,autoStyle,true);
        stylepv={colorPropName c stylePropName ls 'Marker_I',m};
    end
    if multiColor
        stylepv = [stylepv,{'Color', getColor(color, k)}]; %#ok<AGROW> 
    end
    
    h(1,k) = matlab.graphics.chart.primitive.Line('Parent',parax, ...
        ydata{:},xdata{:},...
        stylepv{:},...
        pvpairs{:},extrapairs{k,:});
    h(1,k).assignSeriesIndex();
end
% Create the Surface objects
for k = 1:n
    % detect nan
    nanpos = find(any(isnan(cat(2, getColumn(x,k,xIsNumeric),getColumn(xneg,k,xIsNumeric),getColumn(xpos,k,xIsNumeric),...
                                   getColumn(y,k,yIsNumeric),getColumn(yneg,k,yIsNumeric),getColumn(ypos,k,yIsNumeric)...
                               )),2));
    selpos = {};
    if ~isempty(nanpos)
        selpos{1} = 1:(nanpos(1)-1);
        for p = 2:length(nanpos)
        selpos{p} = (nanpos(p-1)+1):(nanpos(p)-1); %#ok<AGROW> 
        end
        selpos{p+1} = (nanpos(p)+1):sz(1); %#ok<AGROW> 
    end

    % extract data from vectorizing over columns
    if isempty(yneg)  % horizontal
        xdata = {'XData', setPatchVtx(getColumn(x-xneg,k,xIsNumeric), getColumn(x+xpos,k,xIsNumeric), selpos)};
        ydata = {'YData', setPatchVtx(getColumn(y,k,yIsNumeric), getColumn(y,k,yIsNumeric), selpos)};
    else
        xdata = {'XData', setPatchVtx(getColumn(x,k,xIsNumeric), getColumn(x,k,xIsNumeric), selpos)};
        ydata = {'YData', setPatchVtx(getColumn(y-yneg,k,yIsNumeric), getColumn(y+ypos,k,yIsNumeric), selpos)};
    end
    
    stylepv={'FaceColor', h(1,k).Color,'FaceAlpha',0.3, 'EdgeColor','none'};

    h(2,k) = matlab.graphics.primitive.Patch('Parent',parax, ...
        ydata{:},xdata{:},...
        stylepv{:},...
        extrapairs{k,:});
end
% Move plot in front of shading
try %#ok<TRYNC> 
   uistack(h(1,:),'top');
end

if ~hold_state
    set(cax,'Box','on');
end

if nargout>0, hh = h; end

end

%-------------------------------------------------------------------------%
function [pvpairs,args,nargs,msg,orientation,color] = parseargs(args)
% separate pv-pairs from opening arguments
[args,pvpairs] = parseparams(args);

% check existence of color argument
if ~isempty(args)
[color, invalidColors] = matlab.graphics.internal.convertToRGB(args{end});
  if isempty(invalidColors)&&~isempty(color)
    args(end) = [];
  else
    color = [];
  end
end

% Check for LineSpec or Orientation strings
% Allow the orientation flag to occur either before or after the LineSpec
% Allow LineSpec and Orientation to occur at most once each.
validOrientations = {'horizontal','vertical','both'};
orientation = '';
keepArg = true(1,numel(pvpairs));
extraPairs = {};
for a = 1:min(2,numel(pvpairs))
    if isempty(color)
        [color, invalidColors] = matlab.graphics.internal.convertToRGB(pvpairs{a});
        if ~isempty(invalidColors) && ischar(pvpairs{a}) && length(pvpairs{a})>1
            [color, invalidColors] = matlab.graphics.internal.convertToRGB(num2cell(pvpairs{a}));
        end
        if isempty(invalidColors)
            keepArg(a) = false;
        else
            color = [];
        end
    end
    if keepArg(a)
        if matlab.graphics.internal.isCharOrString(pvpairs{a})
            % Check for partial matching of the orientation flag using a
            % minimum of 3 characters.
            tf = strncmpi(pvpairs{a},validOrientations,max(3,numel(pvpairs{a})));
            if isempty(orientation) && any(tf)
                orientation = validOrientations{tf};
                keepArg(a) = false;
            else
                % Check for LineSpec string
                [l,c,m,tmsg]=colstyle(pvpairs{a},'plot');
                if isempty(tmsg) && isempty(extraPairs)
                    keepArg(a) = false;
                    if ~isempty(l)
                        extraPairs = {'LineStyle',l};
                    end
                    if ~isempty(c)
                        extraPairs = [{'Color',c},extraPairs]; %#ok<AGROW>
                    end
                    if ~isempty(m)
                        extraPairs = [{'Marker',m},extraPairs]; %#ok<AGROW>
                    end
                else
                    break;
                end
            end
        else
            % Not a string, so stop looking.
            break
        end
    end
end
if ~isempty(color) && size(color,1)==1
    extraPairs = [extraPairs,{'Color',c}];
    color      = [];
end


linestyleerror = numel(pvpairs)==1;
pvpairs = [extraPairs, pvpairs(keepArg)];
msg = matlab.graphics.chart.internal.checkpvpairs(pvpairs,linestyleerror);
nargs = numel(args);

end

%-------------------------------------------------------------------------%
function val = checkSingleInput(val, sz, argumentName, data)

if isvector(val)
    val = val(:);
end

if ~isempty(sz) && ~isempty(val) && ~isequal(sz, size(val))
    if strcmp(argumentName,'x')
        throwAsCaller(MException(message('MATLAB:errorbar:XDataSizeMismatch')));
    else
        throwAsCaller(MException(message('MATLAB:errorbar:DeltaArgumentSizeMismatch', argumentName)));
    end
end

% Make sure the delta matches the datatype of the data.
if nargin == 4
    if isempty(val)
        if any(strcmp(argumentName,{'err','pos','neg'}))
            % Err is empty, then provide zero-filled matrix.
            val    = data;
            val(:) = 0;
        else
            % Delta is empty, so just make it match the data, regardless of the
            % data provided by the user.
            val = defaultEmptyDelta(data);
        end
    elseif iscategorical(data)
        % Deltas are not supported when data is categorical.
        throwAsCaller(MException(message('MATLAB:errorbar:CategoricalDeltaArgumentNotSupported', argumentName)));
    elseif (isa(data, 'datetime') || isa(data, 'duration')) && ~isa(val, 'duration')
        throwAsCaller(MException(message('MATLAB:errorbar:DeltaArgumentTypeMustBeDuration', argumentName)));
    elseif isnumeric(data) && ~isnumeric(val)
        throwAsCaller(MException(message('MATLAB:errorbar:DeltaArgumentTypeMustBeNumeric', argumentName)));
    end
elseif isempty(val) && strcmp(argumentName,'x')
    % X is empty, then provide sequential matrix.
    val = 1:sz(1);
    val = repmat(val(:),1,sz(2));
end

end

%-------------------------------------------------------------------------%
function delta = defaultEmptyDelta(data)

% Create empty delta based on the data type.
if isa(data, 'datetime') || isa(data, 'duration')
    % duration delta for datetime or duration data
    delta = duration.empty;
else
    % numeric delta for numeric or categorical data
    delta = [];
end

end

%-------------------------------------------------------------------------%
function col = getColumn(val, k, checkData)
    if isempty(val)
        col = val;
    elseif checkData
        col = matlab.graphics.chart.internal.datachk(val(:,k));
    else
        col = val(:,k);
    end
end

%-------------------------------------------------------------------------%
function vtx = setPatchVtx(d1, d2, grp)
    if isempty(grp)
        grp = {1:length(d1)};
    end
    grp(cellfun(@isempty,grp)) = [];

    d1 = d1(:); d2 = d2(:);
    vtx = {};
    for p = 1:length(grp)
        vtx{p} = cat(1, d1(grp{p}), flip(d2(grp{p})) ); %#ok<AGROW> 
    end
    maxNvtx = max(cellfun(@length,vtx),[],'omitnan');
    vtx = cellfun(@(d) cat(1,d,repmat(d(1),maxNvtx-length(d),1)),vtx,'UniformOutput',false);
    vtx = cat(2, vtx{:});
end

%-------------------------------------------------------------------------%
function color = getColor(color, k)
    if size(color,1)>1
        k = mod(k-1,size(color,1))+1;
        color = color(k,:);
    end
end

%-------------------------------------------------------------------------%
function y = issortedarray(x, dim, dir)
    if ~exist('dir','var') || isempty(dir)
        dir = 'ascend';
    end
    if ~exist('dim','var') || isempty(dim)
        dim = 1;
    end
    if isvector(x)
        y = issorted(x,dir);
    else
        y = cellfun(@(r) issorted(r,dir), num2cell(x,dim));
    end
end
