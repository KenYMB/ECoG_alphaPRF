function data = nearlyeq(x,d,mode,k)
% 
% data = NEARLYEQ(X,d)
% data = NEARLYEQ(X,d,mode,n)
%   returns the indices of X nearest to d
% 
% Inputs
%   X     : scalar-array
%   d     : scalar-value or scalar-array
%   n     : returns the first n indices (default = 1)
%           if 'all', it returns all indices where the value is nearest to d
%   mode  :
%       'normal' returns indices of the nearest values (default)
%       'small'  returns indices of the nearest value which is smaller than d
%       'large'  returns indices of the nearest value which is larger than d
%       'both'   returns 'small' in data{:,:,1} and 'large' in data{:,:,2}
% 
% Output
%   data  : scalar-array of indices in X
%           or cell-array containing the indices when 'n~=1' or 'mode=both'
% 

% 2014/07/28 Yuasa
% 2015/04/08 Yuasa - enable scalar-array for d
% 2021/08/26 Yuasa - returns empty when input is nan

narginchk(2,4);

if ~exist('k','var') || isempty(k),    k = 1;   end

if ~exist('mode','var') || isempty(mode),    mode = 'normal'; end

if ~strcmp(mode,'normal') && ~strcmp(mode,'small') && ~strcmp(mode,'large') && ~strcmp(mode,'both'),
    error('nearly equal direction must be ''normal'', ''small'', ''large'' or ''both''.');
end

inpsiz  = size(d);
assert(ndims(d)<=2, 'd can be a scalar or a scalar-array matrix');

if strcmp(mode,'both'),  data = cell([inpsiz 2]);
else                     data = cell(inpsiz);
end

% main
for icol = 1:inpsiz(2)
    for irow = 1:inpsiz(1)

        xs = sort(x(x<=d(irow,icol)),'descend');
        xl = sort(x(x>=d(irow,icol)),'ascend');

        switch mode
            case 'normal'
                if isempty(xs) && isempty(xl),  xs = nan; xl = nan;
                elseif isempty(xs), xs=xl(end);
                elseif isempty(xl), xl=xs(end);
                end
                nearestN = abs(xs(1)-d(irow,icol)) <= abs(xl(1)-d(irow,icol));
                data{irow,icol} = find(x == (xs(1)*nearestN + xl(1)*~nearestN));
            case 'small'
                if min(x(:)) > d(irow,icol) || isempty(xs)
                    data{irow,icol}=zeros(0,1);
                else
                    data{irow,icol} = find(x == xs(1));
                end
            case 'large'
                if max(x(:)) < d(irow,icol) || isempty(xl)
                    data{irow,icol}=zeros(0,1);
                else
                    data{irow,icol} = find(x == xl(1));
                end
            case 'both'
                if min(x(:)) > d(irow,icol) || isempty(xs)
                    data{irow,icol,1}=zeros(0,1);
                else
                    data{irow,icol,1} = find(x == xs(1));
                end
                if max(x(:)) < d(irow,icol) || isempty(xl)
                    data{irow,icol,2}=zeros(0,1);
                else
                    data{irow,icol,2} = find(x == xl(1));
                end
        end

        if ~strcmp(k,'all'),
            if k < 1, k = 1; end
            k = floor(k(1));
            if strcmp(mode,'both'),
                if length(data{irow,icol,1}) > k,   data{irow,icol,1}=data{irow,icol,1}(1:k); end
                if length(data{irow,icol,2}) > k,   data{irow,icol,2}=data{irow,icol,2}(1:k); end
            elseif length(data{irow,icol}) > k
                data{irow,icol}=data{irow,icol}(1:k);
            end
        end
    end
end

if (k(1) == 1) && ~strcmp(mode,'both')
    data = cell2mat(data);
end