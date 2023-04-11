function subjectList = SetSubjectsList(subjectList_fname, varargin)
% SETSUBJECTSLIST returns sujects ID list.
% 
% Usage:
%   subjectList = SetSubjectsList(['subjectlist.tsv'])
%   subjectList = SetSubjectsList(subjectList_fname [,'all'])
%   subjectList = SetSubjectsList(subjectList_fname, select_index)
%   subjectList = SetSubjectsList(subjectList_fname, 'hemi', 'L')
% 
% See also, loadSbjInfo

% 20220224 Yuasa

narginchk(0,inf)
if nargin==1 && isnumeric(subjectList_fname)
    selsbj = subjectList_fname;
    subjectList_fname = '';
elseif nargin==2
    selsbj = varargin{1};
elseif nargin>2
    assert(mod(nargin,2), 'Invalid inputs.');
    selsbj = 'option';
end
SetDefault('subjectList_fname','subjectlist.tsv');
SetDefault('selsbj','all');

%-- Dataset specs
sbjectInfo = loadSbjInfo(subjectList_fname);
if ischar(selsbj)
    switch selsbj
        case 'all'
            selsbj = 1:height(sbjectInfo);
        case 'option'
            selsbj = true(height(sbjectInfo),1);
            nopt   = round((nargin-1)./2);
            for iopt = 1:nopt
                selsbj = selsbj & ismember(sbjectInfo.(varargin{2*iopt-1}),varargin{2*iopt});
            end
        otherwise
            error('Invalid input.');
    end
end
subjectList = sbjectInfo.participant_id(selsbj);
