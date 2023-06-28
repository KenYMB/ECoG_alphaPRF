function    [modelts, datats, cod] = ecog_computePRFtimeseries(stimulus,data,params,opt)
% [modelts, datats, cod] = ecog_computePRFtimeseries(stimulus, data, params, [options])
% [modelts, datats, cod] = ecog_computePRFtimeseries(stimulus, data, result)
%   COMPUTEPRFTIMESERIES rebuild model time series from pRF parameters
% 
% Inputs: 
% - stimulus        = 1xM cell-array of stimulus data.
%                     or CHx1 cell-array of 1xM cell-array of stimulus data.
% - data            = 1xM cell-array of time-series data.
%                     or KxM cell-matrix with bootstrapping.
% - params          = pRF parameter estimates obtained from analyzePRF.
%                     1xCHx7 numeric-matrix, XxCHx7 or XxCHx7xL.
% - options         = a structure with the following fields: 
%       (it's easy to reuse result.options obtained from analyzePRF)
%   - hrf
%   - maxpolydeg
%   - gaussianmode
%   - xvalmode
% - result          = a structure including params and options above.
% 
% Outputs: 
% - modelts         = CH x X x L cell-matrix of model time series estimates.
% - datats          = CH x X x K cell-matrix of data time series.
% - cod             = CH x K x L matrix of coefficient of determination (R^2)
% 
% See also: reconPRFdog(ecog_utils)

% Dependency: <analyzePRF>, SetDefault

% 20200922 Yuasa
% 20200927 Yuasa - consider cross-validation
% 20201027 Yuasa - bug fix
% 20211122 Yuasa - update 'og' model

%% Prepare parameters
%-- Check inputs
narginchk(3,4);
if ~iscell(stimulus),   stimulus = {stimulus};  end
if ~iscell(data),       data = {data};          end
if isstruct(params)
    if isfield(params,'options'),   opt2 = params.options;	end
    params = params.params;
end

%-- Default parametes
SetDefault('opt2.hrf',1,true);
SetDefault('opt.hrf',opt2.hrf,true);
SetDefault('opt2.maxpolydeg',0,true);
SetDefault('opt.maxpolydeg',opt2.maxpolydeg,true);
SetDefault('opt2.gaussianmode','og',false);
SetDefault('opt.gaussianmode',opt2.gaussianmode,false);
SetDefault('opt2.xvalmode',0,false);
SetDefault('opt.xvalmode',opt2.xvalmode,false);

SetDefault('opt2.tr',1,false);
SetDefault('opt.tr',opt2.tr,false);

if opt.xvalmode > 0 && size(params,1)==1
    warning('pRF parameters do not seem to be cross-validated');
    opt.xvalmode = 0;
end

%-- Setup parameters
hrf          = opt.hrf;
maxpolydeg   = opt.maxpolydeg;
gaussianmode = lower(opt.gaussianmode);
tr           = opt.tr;

is3d = size(data{1},4) > 1;
if is3d
  dimdata = 3;
  dimtime = 4;
  xyzsize = sizefull(data{1},3);
else
  dimdata = 1;
  dimtime = 2;
  xyzsize = size(data{1},1);
end
numvxs = prod(xyzsize);

%-- Duplicate stimulus
if ~iscell(stimulus{1}),   stimulus = {stimulus};  end
if length(stimulus)==1,    stimulus = repmat(stimulus,numvxs,1);    end

%-- Check data
res = sizefull(stimulus{1}{1},2);
resmx = max(res);
[nboot1,numruns] = size(data);
nboot2 = size(params,4);
assert(numel(stimulus)==numvxs, '''stimulus'' is unknown format');
assert(length(stimulus{1})==numruns, '''stimulus'' and ''data'' must have same length');

%-- Data length
ntime = zeros(numvxs,numruns);
for pp=1:numruns
    ntime(:,pp) = sum(~isnan(data{1,pp}),dimtime);
end

%-- Setup for cross-validation
switch opt.xvalmode
  case 0
    wantresampleruns = [];
    resampling = {0};
  case 1
    wantresampleruns = 1;
    half1 = copymatrix(zeros(1,numruns),1:round(numruns/2),1);
    half2 = ~half1;
    resampling = {[(1)*half1 + (-1)*half2;
                  (-1)*half1 + (1)*half2]};
    resampling = repmat(resampling,numvxs,1);
  case 2
    wantresampleruns = 0;
    resampling = cell(numvxs,1);
    for vxs = 1:numvxs
    for pp=1:numruns
      half1 = copymatrix(zeros(1,ntime(vxs,pp)),1:round(ntime(vxs,pp)/2),1);
      half2 = ~half1;
      resampling{vxs} = cat(2,resampling{vxs},[(1)*half1 + (-1)*half2;
                                               (-1)*half1 + (1)*half2]);
    end
    end
end
numxval = size(resampling{1},1);
  
%% Compatible
%-- Set typical hrf for empty input
if isempty(hrf)
  hrf = getcanonicalhrf(tr,tr)';
end

%-- Set typical polynomials for empty input
if isempty(maxpolydeg)
  maxpolydeg = cellfun(@(x) round(size(x,dimtime)*tr/60/2),data(1,:));
end
if isscalar(maxpolydeg)
  maxpolydeg = repmat(maxpolydeg,[1 numruns]);
end

%% Prepare for model
%-- Pre-compute cache for faster execution
[~,xx,yy] = makegaussian2d(resmx,2,2,2,2);

%-- Prepare the stimuli for use in the model
stimulusPP = repmat({cell(1,numruns)},numvxs,1);
for vxs=1:numvxs
for pp=1:numruns
  stimulusPP{vxs}{pp} = squish(stimulus{vxs}{pp},2)';  % this flattens the image so that the dimensionality is now frames x pixels
  stimulusPP{vxs}{pp} = [stimulusPP{vxs}{pp} pp*ones(size(stimulusPP{vxs}{pp},1),1)];  % this adds a dummy column to indicate run breaks
end
end

%-- Set model
switch gaussianmode
    case {'dog','gs'}
        modelfun = @(pp,dd) conv2run(modeldogcss(pp(1:5),pp(6:end),dd,res,xx,yy,0,0),hrf,dd(:,prod(res)+1));
    case {'og'}
        modelfun = @(pp,dd) conv2run(posrect(pp(4)) * amppow((dd*[vflatten(placematrix(zeros(res),makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / (2*pi*abs(pp(3))^2))); 0]), posrect(pp(5))),hrf,dd(:,prod(res)+1));
end

%-- Construct projection matrices
polyregressors = repmat({cell(1,numruns)},numvxs,1);
for vxs = 1:numvxs
iboot1=1;
for pp=1:numruns
  if isnan(maxpolydeg(pp))
    polyregressors{vxs}{iboot1,pp} = zeros(ntime(vxs,pp),0);
  else
    polyregressors{vxs}{iboot1,pp} = constructpolynomialmatrix(ntime(vxs,pp),0:maxpolydeg(pp));
  end
end
temp = blkdiag(polyregressors{vxs}{iboot1,:});
cnt = 0;
for pp=1:numruns
    polyregressors{vxs}{iboot1,pp} = temp(cnt+(1:size(polyregressors{vxs}{iboot1,pp},1)),:);
    cnt = cnt + size(polyregressors{vxs}{iboot1,pp},1);
end
end

%-- Figure out test for cross-validation
if ~resampling{1}
    %   testfun =  {@(x) []};
    testfun =  {@(x) catcell(1,x)};
    testfun = repmat(testfun,numvxs,1);
else
    testfun =  cell(numvxs,numxval);
    for vxs = 1:numvxs
        for rr=1:numxval
            testix =  find(resampling{vxs}(rr,:) == -1);
            if wantresampleruns
                testfun{vxs,rr} =  @(x) catcell(1,x(testix));
            else
                testfun{vxs,rr}  = @(x) subscript(catcell(1,x),{testix ':' ':' ':' ':' ':'});
            end
        end
    end
end

%% Loop for voxels (channels)
datats  = cell(numvxs,numxval,nboot1);
modelts = cell(numvxs,numxval,nboot2);
cod     = zeros(numvxs,nboot1,nboot2);
for vxs = 1:numvxs
    tempdata = cellfun(@(x) x(vxs,:)',data,'UniformOutput',0);
    %-- Time series
    for rr=1:numxval
      iboot1 = 1;
      polymatrix = projectionmatrix(feval(testfun{vxs,rr},polyregressors{vxs}(iboot1,:)));
      teststim =  feval(testfun{vxs,rr},stimulusPP{vxs});
      nanplace = isnan(feval(testfun{vxs,rr},tempdata(iboot1,:)));
      for iboot1=1:nboot1
        testdata =  feval(testfun{vxs,rr},tempdata(iboot1,:));
        testdata = testdata(~nanplace);
        datats{vxs,rr,iboot1}  = nan(size(nanplace));
        datats{vxs,rr,iboot1}(~nanplace)  = polymatrix*testdata;
      end
      for iboot2=1:nboot2
        modelpred = feval(modelfun,params(rr,:,vxs,iboot2),teststim(~nanplace,:));
        modelts{vxs,rr,iboot2}  = nan(size(nanplace));
        modelts{vxs,rr,iboot2}(~nanplace)  = nanreplace(polymatrix*modelpred,0,2);
      end
    end

    %-- Variance explained
    for iboot1=1:nboot1
      for iboot2=1:nboot2
        cod(vxs,iboot1,iboot2) = calccod(catcell(1,modelts(vxs,:,iboot2)),catcell(1,datats(vxs,:,iboot1)));
      end
    end
end
          