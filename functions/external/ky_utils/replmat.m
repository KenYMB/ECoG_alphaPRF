function basemat = replmat(basemat,datamat,pos)

% REPLMAT(A,B,position)
%   replace a part of matrix A with matrix B.

% Dependency: modmod

% 20201006 Yuasa

%-- setup
narginchk(2,3);
if ~exist('pos','var') || isempty(pos)
  pos = 1;
end
dimbase = ndims(basemat);
dimdata = ndims(datamat);
sizbase = size(basemat);
sizdata = size(datamat);

assert(dimbase>=dimdata, 'the base matrix must have the same or larger dimensions');
assert(isvector(pos), 'position must be a vector');

%-- set position
if length(pos)==1
    posidx = pos;
    for ii=1:dimbase
        pos(ii) = modmod(posidx,sizbase(ii));
        posidx  = ceil(posidx/sizbase(ii));
    end
end
pos = [reshape(pos,1,[]), ones(1,dimbase-length(pos))];
sizdata = [sizdata, ones(1,dimbase-dimdata)];

%-- figure out indices
idxbase = cell(1,dimbase); 
idxdata = cell(1,dimbase); 
for ii=1:dimbase
    stpos = max(1,pos(ii));
    enpos = min(sizbase(ii),pos(ii)+sizdata(ii)-1);
    idxbase{ii} = stpos:enpos;
    stpos = max(1,-(pos(ii)-2));
    enpos = min(-pos(ii)+sizbase(ii)+1,sizdata(ii));
    idxdata{ii} = stpos:enpos;
end

%-- apply
basemat(idxbase{:}) = datamat(idxdata{:});

