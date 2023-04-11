% ecog_APRFF_INITsub_loadoptions
% General script to load model time-series and pRF data
%   Need to specify parameters listed in "Parameters"

% 20210914
% 20211014 Yuasa: minor update
% 20220202 Yuasa: update for the case alphaType and load options are not consistent
% 20220222 Yuasa: update to avoid using of patient_id

%% Parameters
% Inputs: 
%   subjectList (required if allowmixbeta is true and mixbetaind is not input)
% 
%   alphaType      = 'FaCLb';
%   broadbandType  = 'bbS';
% 
% %   allowbeta    = false;     % if true, load FaCLbB
% %   allowwide    = false;     % if true, load FaCLbW, FaCLbBW
% %     OR
% %   allowmixbeta = false;     % if true, set beta options automatically
% %     mixbetaind              % you can manually set beta options
% %   
% %   allowlag     = false;     % if true, load 'lag' data
% %   linalpha     = false;     % if true, load FaCRb
% %   noACorrect   = false;     % if true, load FaLb

%-- check empty
AutoFILLalphaType  = ~exist('alphaType','var')||isempty(alphaType);
UPDATEcomputeOpt  = false;
hasmixbetaind = exist('mixbetaind','var') && ~isempty(mixbetaind);

%-- set default values (if possible, estimate options)
if ~AutoFILLalphaType
    UPDATEcomputeOpt = UPDATEcomputeOpt | UPDATEcomputeOptFnc('allowmixbeta', iscell(alphaType) && length(alphaType)>1 && ~all(matches(alphaType,alphaType{1})),1);
    UPDATEcomputeOpt = UPDATEcomputeOpt | UPDATEcomputeOptFnc('allowwide', any(endsWith(alphaType,{'W','Wlag'})));
    UPDATEcomputeOpt = UPDATEcomputeOpt | UPDATEcomputeOptFnc('allowbeta', any(endsWith(alphaType,{'B','BW','Blag','BWlag'})));
    UPDATEcomputeOpt = UPDATEcomputeOpt | UPDATEcomputeOptFnc('noACorrect', any(~startsWith(alphaType,{'aC','FaC'})));
    UPDATEcomputeOpt = UPDATEcomputeOpt | UPDATEcomputeOptFnc('linalpha', any(startsWith(alphaType,{'aR','FaR','aCR','FaCR'})));
    SetDefault('allowlag',false);
    if UPDATEcomputeOpt
        warning('Load options are updated based on input alphaType.');
    end
else
    SetDefault('allowmixbeta',false);
    SetDefault('allowbeta',false);
    SetDefault('allowwide',false);
    SetDefault('allowlag',false);
    SetDefault('linalpha',false);
    SetDefault('noACorrect',false);
    if allowmixbeta
        if ~allowbeta && ~allowwide   % consider both are true when both allowbeta & allowwide are false
            allowbeta  = true;
            allowwide  = true;
        end
        if ~hasmixbetaind
            mixbetaind = alphaFitTypes(subjectList,'index'); % 0=nobeta, 1=beta, 2=wide, 3=betawide
        end
        if all(mixbetaind == mixbetaind(1))  % correct if all subjects use the same options
%             allowmixbeta = false;
            allowbeta    = mod(mixbetaind(1),2)==1;
            allowwide    = fix(mixbetaind(1)./2)>=1;
        end
    end
    if noACorrect
        allowbeta    = false;
        allowwide    = false;
        allowmixbeta = false;
    end
end
    
%%% process options
%-- broadband type
SetDefault('broadbandType','bbS');
%-- basic alpha type
if linalpha && noACorrect
    SetDefault('alphaType','FaRb');
elseif linalpha 
    SetDefault('alphaType','FaCRb');
elseif noACorrect
    SetDefault('alphaType','FaLb');
else
    SetDefault('alphaType','FaCLb');
end
%-- beta fit option in alpha fitting
if AutoFILLalphaType
    if allowmixbeta
        alphaType = repmat({[alphaType]},numel(subjectList),1);
        alphaType(mixbetaind==1) = {[alphaType{find(mixbetaind==1,1)} 'B']};
        alphaType(mixbetaind==2) = {[alphaType{find(mixbetaind==2,1)} 'W']};
        alphaType(mixbetaind==3) = {[alphaType{find(mixbetaind==3,1)} 'BW']};
    else
        if allowbeta, alphaType = [alphaType 'B']; end
        if allowwide, alphaType = [alphaType 'W']; end
    end
end
%-- lag option in regression
if allowlag
    if ~endsWith(broadbandType,'lag'), broadbandType = [broadbandType 'lag']; end
    if ~iscell(alphaType)
      if ~endsWith(alphaType,'lag'),     alphaType = [alphaType 'lag'];         end
    else
      for ii=1:length(alphaType)
        if ~endsWith(alphaType{ii},'lag'), alphaType{ii} = [alphaType{ii} 'lag'];  end
      end
    end
end
%-- feedback to allowlag
if endsWith(broadbandType,'lag') && all(endsWith(alphaType,'lag'),'all')
     allowlag = true;
end

%-- reset mixbetaind
if ~hasmixbetaind
    clear mixbetaind
end

%% R2mode
if exist('usefulltsR2','var')&&exist('usefulltsxR2','var')
    usefllts = usefulltsR2 | usefulltsxR2;
    %-- name label
    if usefllts
        R2mode = '-fullts';
    else
        R2mode = '';
    end
end

%% subfunction
function    isupdated = UPDATEcomputeOptFnc(VARNAME,VALUE,KEEPVAR)
if ~exist('KEEPVAR','var')||isempty(KEEPVAR)
    KEEPVAR = false;
end
prevval   = [];
if evalin('caller',sprintf('exist(''%s'',''var'')',VARNAME))   % check existence
    prevval   = evalin('caller',VARNAME);
end
%-- Update values
if KEEPVAR&&~isempty(prevval)
    VALUE = prevval;
end
assignin('caller',VARNAME,VALUE);
isupdated = ~isempty(prevval) && ~isequaln(VALUE,prevval);
end
