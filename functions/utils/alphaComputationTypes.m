function alphaTypes = alphaComputationTypes(selTStype)

% The list of target band names for alpha.
%   'a'    = spectral power change at alpha frequnecy
%   'aC'   = estimated alpha power change with gaussian fitting
%            (estimateIAF = false, allownegfit = false in ecog_prf_fitalpha)
%   'aCR'  = inverse value of 'aC'
%   'aCL'  = 'aC' in log scale
%   'F*'   = 'with fixing alpha frequency at IAF
%            (estimateIAF = true in ecog_prf_fitalpha)
%   '*b'   =  with allowing negative coefficients
%            (allownegfit = true in ecog_prf_fitalpha)
%   '*B'   =  with allowing beta fitting
%            (allowbetafit = true in ecog_prf_fitalpha)
%   '*BW'  = 'aC' with allowing beta fitting in wide frequency range
%            (allowbetafit = true, allowwidefit = true in ecog_prf_fitalpha)
%   '*lag' = 'aC' on regressed data with allowing 0.01s lags
%            (allowlag = true in ecog_prf_regressData)

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
params  = {{'','C'},{'','R','L'}};
nparams = length(params);
nList   = cell(1,nparams);
[nList{:}] = ndgrid(params{:});
nList      = reshape(cat(nparams+1,nList{:}),[],nparams);
tstype = cell(0);
for ii=1:size(nList,1)
    tstype{1,ii} = cat(2,nList{ii,:});
end

if ~ismember('*',selTStype)
    tstype = tstype(ismember(tstype,selTStype));
end

%-- make list
params  = {{'','F'},{'a'},tstype,{'','b'},{'','B','BW'},{'','lag'}};
nparams = length(params);
nList   = cell(1,nparams);
[nList{:}] = ndgrid(params{:});
nList      = reshape(cat(nparams+1,nList{:}),[],nparams);
alphaTypes = cell(0);
for ii=1:size(nList,1)
    alphaTypes{1,ii} = cat(2,nList{ii,:});
end

end
