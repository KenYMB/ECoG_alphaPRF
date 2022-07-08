function    adjustplot(h)
% ADJUSTPLOT(h)
%   bounds OuterPosition of each axis into [0 1]

% 20210105 Yuasa

if ~exist('h','var')
    h = get(groot,'CurrentFigure');
end

if isgraphics(h)&&ismember(h.Type,{'figure'})
    curax = get(h,'Children');  drawnow;
    for ii = 1:length(curax)
        curax(ii).OuterPosition = min(max(curax(ii).OuterPosition,0),1);
        drawnow;
    end
end
drawnow;