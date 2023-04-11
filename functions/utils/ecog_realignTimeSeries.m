function [modeldata] = ecog_realignTimeSeries(modeldata,repstim)

% Description: 
%
% [modeldata_sbj] = ecog_realignTimeSeries(modeldata_sbj)
% [modeldata_sbj] = ecog_realignTimeSeries(modeldata_sbj,replicate_stimlus=True)
% 
% ECOG_REALIGNTIMESERIES insert NaNs at the end of each time series to
% have the same time length across input data. Useful to merge data time
% series across subjects.
% If replicate_stimlus is true, the 'stimulus' field is replicated to have
% the same length as the numbers of channels.
% 
% 
% Input & Output
% - modeldata_sbj   = Nx1 cell-array of time-series structure for each subject
% 
% See also: ECOG_PRF_CONSTRUCTTIMESERIES, ECOG_REARRANGEPRF

% Dependency: <analyzePRF>, SetDefault, replmat, int2ordinal

% 20201005 Yuasa

%% parameter setting
SetDefault('repstim',false);
%-- convert structure-array to cell-array
if isstruct(modeldata)
    modeldata = mat2cell(modeldata,ones(size(modeldata)));
elseif ~iscell(modeldata)
    modeldata = {modeldata};
end

%-- check data validity & get maximum data size
ndat = numel(modeldata);
nrun = 1;
nitr = 1;
for ii = 1:ndat
    assert(all(ismember({'datats','stimulus'},fieldnames(modeldata{ii}))), ...
        'Invalid input (%s data)',int2ordinal(ii));
    assert(~iscell(modeldata{ii}.stimulus{1}),...
        '%s might be already applied.',mfilename);
    assert(size(modeldata{ii}.datats,2)==size(modeldata{ii}.stimulus,2), ...
        'Invalid input (%s data)',int2ordinal(ii));
    nrun = max(nrun, size(modeldata{ii}.datats,2));
    nitr = max(nitr, size(modeldata{ii}.datats,1));
end
ntim = ones(1,nrun);
for ii = 1:ndat
    for jj = 1:size(modeldata{ii}.datats,2)
        assert(size(modeldata{ii}.datats{1,jj},2)==size(modeldata{ii}.stimulus{1,jj},3), ...
            'Invalid input (%s data, %s run)',int2ordinal(ii),int2ordinal(jj));
        kk=1;
        ntim(kk,jj) = max(ntim(kk,jj), size(modeldata{ii}.datats{kk,jj},2));
    end
end

%% insert NaNs
for ii = 1:ndat
    [citr,ctim] = size(modeldata{ii}.datats);
    sizdat = size(modeldata{ii}.datats{1});
    sizstm = size(modeldata{ii}.stimulus{1});
    modeldata{ii}.datats = placematrix2(cell(nitr,nrun),modeldata{ii}.datats,1);
    modeldata{ii}.stimulus = placematrix2(cell(1,nrun),modeldata{ii}.stimulus,1);
    for jj = 1:nrun
        dummydat = cast(nan(copymatrix(sizdat,2,ntim(jj))),class(modeldata{ii}.datats{1,jj}));
        dummystm = cast(false(copymatrix(sizstm,3,ntim(jj))),class(modeldata{ii}.stimulus{1,jj}));
        for kk = 1:nitr
          if jj<=ctim && kk<=citr
            modeldata{ii}.datats{kk,jj} = replmat(dummydat,modeldata{ii}.datats{kk,jj},1);
          else
            modeldata{ii}.datats{kk,jj} = dummydat;
          end
        end
        if jj<=ctim
          modeldata{ii}.stimulus{1,jj} = replmat(dummystm,modeldata{ii}.stimulus{1,jj},1);
        else
          modeldata{ii}.stimulus{1,jj} = dummystm;
        end
    end
end

%% replicate stimulus
if repstim
    for ii = 1:ndat
        cchn = size(modeldata{ii}.datats{1},1);
        modeldata{ii}.stimulus = repmat({modeldata{ii}.stimulus},cchn,1);
    end
end
