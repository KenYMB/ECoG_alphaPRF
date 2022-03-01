function saveauto(filename, varargin)

% SAVEAUTO is alternative function of SAVE automatically adding '-v7.3'
% option if the data size will be larger than 2GB.
% 
% see also, save

% 20171212 Yuasa
% 20190816 Yuasa: bug fix on '-struct' option

%-- parameter
savethresh = 2^31;

%-- separate datas and options
datanames = varargin;
optnames  = varargin;
isv73opt  = false;
isstruct  = false;
isregexp  = false;
for ilp=(nargin-1):-1:1
    if isempty(varargin{ilp})               % skip empty
        datanames(ilp) = [];
        optnames(ilp)  = [];
    else
        if strcmp(varargin{ilp}(1),'-')         % check if the input is option
            if strcmpi(varargin{ilp},'-v7.3')
                isv73opt   = true;
            elseif strcmpi(varargin{ilp},'-struct')
                assert(ilp<nargin-1,message('MATLAB:save:missingStructArg'));
                isstruct   = true;
                strctname        = varargin{ilp+1}; % get structure name following '-struct'
                datanames(ilp+1) = [];              % delete structure name from datanames
            elseif strcmpi(varargin{ilp},'-regexp')
                isregexp   = true;
            end
            datanames(ilp) = [];                % delete option from datanames
        else
            optnames(ilp)  = [];                % delete variable from optnames
        end
    end
end

%-- main
if isv73opt                                     % passthrough
    evalin('caller',[sprintf('save(''%s''',filename),  sprintf(',''%s''',varargin{:}), sprintf(');')]);
else
    wkspvar    = evalin('caller','who');
    if isempty(datanames)                       % check exist of dataname or filedname input
        datanames = {'.*'};
    elseif ~isregexp                            % convert into regexp
        datanames = regexptranslate('wildcard',datanames);
    end
    
    %%% -struct
    if isstruct
        assert(any(strcmp(wkspvar,strctname)),message('MATLAB:save:structArgNotScalarStruct'));
        %-- check each fileds
        ilp = 1;
        fldlist = evalin('caller', sprintf('fieldnames(%s)',strctname));
        while ~(isv73opt || ilp>length(datanames))
            ilpsub = 1;
            fldlistsub = regexp(fldlist,['^' datanames{ilp} '$'],'match');
            fldlistsub = cat(1,fldlistsub{:});
            while ~(isv73opt || ilpsub>length(fldlistsub))
                tmpfld = evalin('caller',[strctname '.' fldlistsub{ilpsub}]);
                isv73opt = getfield(whos('tmpfld'),'bytes') >= savethresh;
                ilpsub = ilpsub + 1;
            end
            ilp = ilp + 1;
        end
    %%% not -struct
    else
        %-- check each datanames
        ilp = 1;
        while ~(isv73opt || ilp>length(datanames))
            ilpsub = 1;
            tmplist = regexp(wkspvar,['^' datanames{ilp} '$'],'match');
            tmplist = cat(1,tmplist{:});
            while ~(isv73opt || ilpsub>length(tmplist))
                isv73opt = evalin('caller',sprintf('getfield(whos(''%s''),''bytes'')',tmplist{ilpsub})) >= savethresh;
                ilpsub = ilpsub + 1;
            end
            ilp = ilp + 1;
        end
    end
    
    if isv73opt
        evalin('caller',[sprintf('save(''%s''',filename),  sprintf(',''%s''',varargin{:}), sprintf(',''-v7.3'');')]);
    else
        evalin('caller',[sprintf('save(''%s''',filename),  sprintf(',''%s''',varargin{:}), sprintf(');')]);
    end
end

