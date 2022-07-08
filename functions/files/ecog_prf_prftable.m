function [prf_table] = ecog_prf_prftable(channels,varargin)
% [prf_table] = ecog_prf_prftable(channels,prf_params)
% [prf_table] = ecog_prf_prftable(channels,channel_fields,prf_params)
% [prf_table] = ecog_prf_prftable(channels,prf_params1,postfix1,prf_params2,postfix2,...)
%    summarize pRF results in table
% 
% Outputs:
%   prf_table       = table including channel names and pRF results
% 
% Inputs:
%   [pRF results]
%   channels        = channel table
%   prf_params      = structure of pRF parameters
% 
%   [options]
%   channel_fields  = field name which are copied from channels 
%                     'name' and 'subject_name' are automatically copied
%   postfix         = identifier to merge multiple prf_params
% 
% Example:
%   prf_table = ecog_prf_prftable(channels,'wangarea',prf_all_bb,'bb',prf_all_a,'a');

% 20220518 Yuasa

%% interpret inputs
cpfields = {'name', 'subject_name','hemisphere'};
if nargin > 1 && (ischar(varargin{1})||iscellstr(varargin{1})||isstring(varargin{1}))
    cpfields = horzcat(cpfields,...
                        reshape(cellstr(varargin{1}),1,[]));
    varargin(1) = [];
end

prf     = {};
postfix = {};
ll = 1;
while ll <= length(varargin)
    prf(end+1) = varargin(ll);
    if ll < length(varargin) && ischar(varargin{ll+1})
        postfix{end+1} = sprintf('_%s',varargin{ll+1});
        ll = ll + 2;
    else
        postfix{end+1} = '';
        ll = ll + 1;
    end
end

%% prepare table
chanflds = fieldnames(summary(channels));
cpfields = intersect(cpfields,chanflds,'stable');

prf_table = channels(:,cpfields);

%% copy prf results
prfparams = {'xval','ecc','ang','rfsize','expt','gain'};

for prfprm = prfparams
    for ll = 1:length(prf)
        prfflds = fieldnames(prf{ll});
        oupfld  = sprintf('%s%s',prfprm{:},postfix{ll});
        inpfld  = prfprm{:};
        if strcmpi(inpfld,'xval') && ~ismember('xval',prfflds) && ismember('aggregatedtestperformance',prfflds)
            inpfld = 'aggregatedtestperformance';
        end
        if ismember(inpfld,prfflds)
            assert(size(prf{ll}.(inpfld),1)==height(channels),...
                '''prf_params'' must consist of the same channels with ''channels''');
            prf_table.(oupfld) = prf{ll}.(inpfld);
        end
    end
end
