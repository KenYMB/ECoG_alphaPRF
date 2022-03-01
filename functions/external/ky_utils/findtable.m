function I = findtable(T,varargin)

% I = FINDTABLE(T,VariableName,Variables)
%   returns row indices where T.VariableName is found in Variables.
% 
% I = FINDTABLE(T,VariableName1,Variables1,VariableName2,Variables2,...)
%   returns row indices where all VariableNames are found in respective Variables.

% 20210105 Yuasa

%-- check inputs
assert(istable(T),'First argument must be a table');
assert(logical(mod(nargin,2)),'VariableName and Variables must be specified as pairs');

%-- reshape Variables into columns
for vars = 2:2:length(varargin)
    varargin{vars} = reshape(varargin{vars},[],1);
end

%-- make reference table
T2 = table(varargin{2:2:end},'VariableNames',varargin(1:2:end));

%-- remove unnecessary Variables
[~,TV1,TV2] = intersect(T.Properties.VariableNames,T2.Properties.VariableNames);
if length(TV2) == size(T2,2)
    T1 = T(:,TV1);
    
%-- find
    I  = find(ismember(T1,T2));
    
else
%-- output empty if T does not have VariableName
    I = zeros(0,1);
end