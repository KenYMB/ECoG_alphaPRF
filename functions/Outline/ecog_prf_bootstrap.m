function    [out] = ecog_prf_bootstrap(data,N,opts)
% [data] = ecog_prf_bootstrap(data,N,opts)


% 20200407 Yuasa


%% Set options
%--Define inputs 
narginchk(1,3);
SetDefault('N',100);
% <opts>
SetDefault('opts.issave',false);
SetDefault('opts.outputDir',fullfile(analysisRootPath, 'Data'));
SetDefault('opts.exclude',{},'cell');

%% Loop across subjects
cellinput = iscell(data);
if ~cellinput,  data = {data};  end
if iscolumn(data)
    out = cell([length(data),N]);
else
    out = cell([size(data),N]);
end

for ii = 1:numel(data)
    idat = data{ii};
    
    %-- Find fields depending event numbers
    eventsize = height(idat.events);
    flds      = fieldnames(idat);
    flds(structfun(@(x) all(size(x)~=eventsize),idat)) = [];
    flds(ismember(flds,opts.exclude)) = [];
    
    %-- Get stimulus list & prepare for bootstrap
    stmlist     = findgroups(idat.events.trial_name);
    bootevents  = zeros(eventsize,N);
    for jj = reshape(unique(stmlist),1,[])
        eventidx = find(stmlist==jj);
        eventnum = numel(eventidx);
        bootidx  = eventidx(floor(rand(length(eventidx),N).*eventnum)+1);
        
        bootevents(eventidx,:) = bootidx;
    end
    
    %%% Bootstrap
    for bb = 1:N
        bdat = idat;
        bootidx = bootevents(:,bb);
        for jj = reshape(flds,1,[])
            eventdim = sort(find(size(bdat.(jj{:}))==eventsize),'descend');
            smplstr  = repmat(':,',1,ndims(bdat.(jj{:}))); smplstr(end) = [];
            for kk = 2.*eventdim-1
                smplstr = [smplstr(1:(kk-1)) 'bootidx' smplstr((kk+1):end)];
            end
            
            bdat.(jj{:}) = eval(['bdat.(jj{:})(' smplstr ')']);
        end
        
        %-- Output
        out{ii,bb} = bdat;
    end
end
   
if ~cellinput,  out = out{1};  end         
            