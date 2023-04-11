function alphaFitTypes = alphaFitTypes(subjectList,outputmode)

% alphaFitTypes = ALPHAFITTYPES(subjectList, ['index'])
%   outputs how to compute alpha for each subject as the following index,
%   by loading 'alphafittype.tsv'.
%       0: normal
%       1: beta
%       2: wide
%       3: betawide
% 
% alphaFitTypes = ALPHAFITTYPES(subjectList, 'name')
%   outputs the names of how to compute alpha for each subject.
% 

% 20220222 Yuasa

%-- Default Value
narginchk(1,2);
if nargin < 2 || isempty(outputmode)
    outputmode = 'index';
end

%-- Load table
ATypes  =readtable('alphafittype.tsv', 'FileType', 'text');
ATypes.typeid(:)                                     = 0;
ATypes.typeid(ismember(ATypes.recommend,'beta'))     = 1;
ATypes.typeid(ismember(ATypes.recommend,'wide'))     = 2;
ATypes.typeid(ismember(ATypes.recommend,'betawide')) = 3;

%-- Find subjects
if ~iscell(subjectList)&&~isstring(subjectList),  subjectList = {subjectList};  end
[~,idx1,idx2] = intersect(subjectList,ATypes.participant_id,'stable');

%-- Outputs
switch outputmode
    case 'index'
        alphaFitTypes = zeros(size(subjectList));
        alphaFitTypes(idx1) = ATypes.typeid(idx2);
    case 'name'
        alphaFitTypes = repmat({'normal'},size(subjectList));
        alphaFitTypes(idx1) = ATypes.recommend(idx2);
    otherwise
        error('Unknown mode');
end