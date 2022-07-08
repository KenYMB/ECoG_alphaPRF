function    varargout = getelement(datamat, varargin)
% d  = GETELEMENT(A,i1,i2,...)
%   getelement(A,i1,i2) is same as A(i1,i2)
%   getelement(A,'end',':') is same as A(end,:)
%   getelement(A,i1,i2,'cell') is same as A{i1,i2}
%   getelement(A,i1,':','cell','cat') is same as [A{i1,:}]

% 20170825 Yuasa
% 20191227 Yuasa: fix for the case the number of arguments is less than the
%                 dimension of data matrix

if nargin < 2
    vargin{1} = ':';
end

%-- check cell
iscell =  cellfind(varargin,'cell');
if ~isempty(iscell)
    varargin(iscell) = [];
    iscell = true;
else
    iscell = false;
end

%-- check cat
iscat =  cellfind(varargin,'cat');
if ~isempty(iscat)
    varargin(iscat) = [];
    iscat = true;
else
    iscat = false;
end

%-- reshape
dimelem = nargin -1 -iscell -iscat;
dimmat  = ndims(datamat);
if ~isvector(datamat) &&  dimelem < dimmat
    if dimelem == 1
        newsiz = [numel(datamat) 1];
    else
        newsiz = [size(datamat,1:(dimelem-1)) prod(size(datamat,dimelem:dimmat))];
    end
    datamat = reshape(datamat,newsiz);
end

%-- find 'end'
for iinp = 1:length(varargin)
    if ischar(varargin{iinp}) && ~strcmp(varargin{iinp},':')
        if iinp > ndims(datamat),   dimsiz = '1';
        else,                       dimsiz = num2str(size(datamat,iinp));
        end
        varargin{iinp} = strrep(varargin{iinp},'end',dimsiz);
        varargin{iinp} = eval(varargin{iinp});
    end    
end

%-- main
if iscell
  if iscat
    varargout{1} = [datamat{varargin{:}}];
  else
    varargout    = datamat(varargin{:});
  end
else
    varargout{1} = datamat(varargin{:});
end    