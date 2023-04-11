%% ECoG Alpha pRF
% ecog_APRF_01a_preprocessing
%   load ECoG data from BIDS files

% Check if compuatation is completed
% # IDAT=4801
% # for ((i = IDAT ; i <= IDAT+30 ; i++)); do FILE=log/error_fitPRFboot-prep-${i}.txt; if [ -f "$FILE" ] ; then echo $i; tail $FILE; fi; done;
% # for ((i = IDAT ; i <= IDAT+30 ; i++)); do FILE=log/out_fitPRFboot-prep-${i}.txt; if [ -f "$FILE" ] ; then echo $i; tail -n 2 $FILE; fi; done;
% #
% # IDAT=4000
% # for ((j=IDAT; j<=IDAT+900 ; j=j+100)); do for ((i=j; i<=j+30; i++)); do FILE=log/error_fitPRFboot-prep-${i}.txt; if [ -f "$FILE" ] ; then echo $i; tail $FILE; fi; done; done

% 20210907 Yuasa: run for allbeta options
% %% without ERP %%

% See also: ecog_APRFF_01d0_analyzePRFboot_prep

%% Define paths and dataset
checkPath;
%-- Subjects
SetDefault('subjectList_fname','subjectlist.tsv');
SetDefault('selsbj','all');
subjectList = SetSubjectsList(subjectList_fname, selsbj);

%-- modifiable parameters
SetDefault('allowlag',false);
SetDefault('average','runs');
SetDefault('allowbetafit',ismember(alphaFitTypes(subjectList,'name'),{'betawide','beta'}));
SetDefault('allowwidefit',ismember(alphaFitTypes(subjectList,'name'),{'betawide','wide'}));
SetDefault('bootstrapindex',nan);

%% Merge fit alpha
%-- parameters
gammafit       = false;
estimateIAF    = true;
allownegfit    = true;
numbootstrap   = 1000;

%-- file location
outputDir = 'Spectrum';
fileid    = ['freq_spectra-params-boot'];

%-- file name parameters
postfix = '';
%%-- lag
addlag = allowlag && ~contains(fileid,'regresslag');
if addlag,   postfix = sprintf('%s_regresslag',postfix);	end
%%-- parameters 
fitparams = '';
addbeta = allowbetafit && ~contains(fileid,{'_beta'});
addwide = allowwidefit && ~contains(fileid,{'wide_'}) && ~endsWith(fileid,{'wide'}) ;
if addbeta,         fitparams = sprintf('%sbeta',fitparams); end
if addwide,         fitparams = sprintf('%swide',fitparams); end
if ~estimateIAF, 	fitparams = sprintf('%s-freeAlpha',fitparams); end
if ~allownegfit,  	fitparams = sprintf('%s-oneside',fitparams); end
if ~gammafit,       fitparams = sprintf('%s-nogammafit',fitparams); end
fitparams    = regexprep(fitparams,'^-*','');
postfix      = sprintf('%s_%s_avg-%s',postfix,fitparams,average);

%-- Run for subjects
outputPth = SetDefaultAnalysisPath('DAT',outputDir);
for ii = 1:length(subjectList)
    %-- file name
    subject      = subjectList{ii};
    filenamebase = sprintf('%s_%s-*%s.mat', subject,fileid,postfix);
    filenames    = dir(fullfile(outputPth, filenamebase));
        filenames       = {filenames(:).name}';
        fileinvalid     = cellfun(@isempty,regexp(filenames,[fileid '-\d+' postfix],'match','once'));
        filenames(fileinvalid) = [];
        filenums        = str2double(regexp(regexp(filenames,[fileid(end-5:end) '-\d+' postfix(1:5)],'match','once'),'\d+','match','once'));
        [~, fileorder]  = sort(filenums(~isnan(filenums)));
        filenames       = filenames(fileorder);
    filename     = sprintf('%s_%s%s.mat', subject,fileid,postfix);
        filepaths    = fullfile(outputPth,filenames);
        filepath     = fullfile(outputPth,filename);
        
    %-- check availability
    if isempty(filenames)
        warning('No separate file is found. Skip merging...');
    else
        %-- load files & merge
        for jj = 1:length(filepaths)
            fprintf(1,'Loading %s from %s\n',filenames{jj},outputPth);
            tmp = load(filepaths{jj});
            if jj==1
                iparams = tmp;
            else
                iparams.resamp_parms = ...
                    cellfun(@(B,C) cat(3,B,C), iparams.resamp_parms, tmp.resamp_parms, 'UniformOutput', false);
            end
        end
        
        %-- check data size
        bootnums = cellfun(@(C) size(C,3),iparams.resamp_parms);
        dltflg   = false;
        if max(bootnums) ~= min(bootnums)
            warning('Some electrodes missed bootstrap data: Max# %d, Min# %d.', max(bootnums),min(bootnums));
        elseif min(bootnums) < numbootstrap
            warning('Total number of bootstrapping (%d) is less than expected (%d).', min(bootnums), numbootstrap);
        else
            if max(bootnums) > numbootstrap
                warning('Total number of bootstrapping (%d) is more than expected (%d).', max(bootnums), numbootstrap);
                bootselect = randi(max(bootnums),numbootstrap,1);
                iparams.resamp_parms = cellfun(@(P) P(:,:,bootselect),iparams.resamp_parms,'UniformOutput',false);
                warning('Deleting excessive iterations.');
            end
            %-- delete separated data
            dltflg = true;
        end
        
        %-- save merged data
        fprintf(1,'Loaded data is merged into %s\n',filename);
        if exist(filepath,'file')
            warning('%s is overwitten.',filepath);
        end
        saveauto(filepath,'-struct','iparams');
        
        %-- delete separated data
        if dltflg
            delete(filepaths{:});
            warning('Individual data files are deleted');
        end
        
    end
end
