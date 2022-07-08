function    realignYLabels(h)
% REALIGNYLABELS(h)
%   realigns horizontal positions of y labels to the left end label.
%   h is an array of axes or text.

% 20220606 Yuasa

if ~isempty(h)
    %-- Check input
    if iscell(h),  h = vertcat(h{:});   end
    assert(all(ishandle(h),'all'),message('MATLAB:class:InvalidHandle'));
    if isa(h,'matlab.graphics.axis.Axes') && ismember('YLabel',fieldnames(h))
        h = vertcat(h.YLabel);
    end
    
    %-- Get position of the left end label
    pause(0.1);     % need to refresh axes
    ylpos  = vertcat(h.Position);
    ylxlim = arrayfun(@(x) x.Parent.XLim,h,'UniformOutput',false);
    ylxlim = vertcat(ylxlim{:});
    ylrelx = (ylxlim(:,1)-ylpos(:,1))./diff(ylxlim,[],2);
    [maxrelx, whichyl] = max(ylrelx);
    ylxloc = ylxlim(:,1) - diff(ylxlim,[],2) .* maxrelx;
    ylpos(:,1) = ylxloc;
    
    %-- Apply realign
    yy =1:length(h);  yy(whichyl) = [];
    for yy=yy
        h(yy).Position = ylpos(yy,:);
    end
end