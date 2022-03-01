function broadbandTypes = broadbandComputationTypes(selTStype)

% The list of target band names for broadband.
%   'bb'      = broadband estimated with linear regression
%               (gammafit = true in ecog_prf_fitalpha)
%   'bbS'     = broadband estimated as mean power
%               (gammafit = false in ecog_prf_fitalpha)
%   'bbL'     = broadband estimated in low frequency range
%   '*lag'  = 'bbS' on regressed data with allowing 0.01s lags
%               (allowlag = true in ecog_prf_regressData)

% 20211117 Yuasa

%-- set input
SetDefault('selTStype','*',true,'cell');
if ~iscellstr(selTStype)
    if numel(selTStype)==1 && isempty(selTStype{1})
        selTStype = {''};
    else
        error('Invalid input');
    end
end

%-- set selection
tstype = {'','S','L'};

if ~ismember('*',selTStype)
    tstype = tstype(ismember(tstype,selTStype));
end

%-- make list
params  = {{'bb'},tstype,{'','lag'}};
nparams = length(params);
nList   = cell(1,nparams);
[nList{:}] = ndgrid(params{:});
nList      = reshape(cat(nparams+1,nList{:}),[],nparams);
broadbandTypes = cell(0);
for ii=1:size(nList,1)
    broadbandTypes{1,ii} = cat(2,nList{ii,:});
end

end