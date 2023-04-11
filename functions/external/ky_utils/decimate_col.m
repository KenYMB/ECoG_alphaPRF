function odata = decimate_col(idata,r,nfilt,option,dim,nanflag)
%DECIMATE_COL Resample a data vector at a lower rate after lowpass filtering.
%   For matrices or N-D arrays, decimate the elements along the first array
%   dimension.
%   Y = DECIMATE_COL(X,R) resamples the sequence in matrix X at 1/R times the
%   original sample rate.  The resulting resampled matrix Y is R times
%   shorter, i.e., SIZE(Y,1) = CEIL(SIZE(X,1)/R). By default, DECIMATE_COL
%   filters the data with an 8th order Chebyshev Type I lowpass filter with
%   cutoff frequency .8*(Fs/2)/R, before resampling.
%
%   Y = DECIMATE_COL(X,R,N) uses an N'th order Chebyshev filter.  For N greater
%   than 13, DECIMATE will produce a warning regarding the unreliability of
%   the results.  See NOTE below.
%
%   Y = DECIMATE_COL(X,R,'FIR') uses a 30th order FIR filter generated by
%   FIR1(30,1/R) to filter the data.
%
%   Y = DECIMATE_COL(X,R,N,'FIR') uses an Nth FIR filter.
% 
%   Y = DECIMATE_COL(X,R,[],dim)
%   Y = DECIMATE_COL(X,R,N,dim)
%   Y = DECIMATE_COL(X,R,'FIR',dim)
%   Y = DECIMATE_COL(X,R,N,'FIR',dim)
%       decimate X along dimension dim. 
% 
%   Y = DECIMATE_COL(...,'omitnan')
%       ignores all NaN values in the calculation
%
%   Note: For better results when R is large (i.e., R > 13), it is
%   recommended to break R up into its factors and calling DECIMATE several
%   times.
%
%   EXAMPLE: Decimate a signal by a factor of four
%   t = 0:.00025:1;  % Time vector
%   x = [sin(2*pi*30*t) + sin(2*pi*60*t);
%        sin(2*pi*40*t) + 0.4*sin(2*pi*80*t)]';
%   y = decimate_col(x,4);
%   subplot(1,2,1);
%   stem(x(1:120,:)), axis([0 120 -2 2])   % Original signal
%   title('Original Signal')
%   subplot(1,2,2);
%   stem(y(1:30,:))                        % Decimated signal
%   title('Decimated Signal')
%
%   See also DECIMATE, DOWNSAMPLE, INTERP, RESAMPLE, FILTFILT, FIR1, CHEBY1.

%   Author(s): L. Shure, 6-4-87
%              L. Shure, 12-11-88, revised
%   Copyright 1988-2018 The MathWorks, Inc.

%   References:
%    [1] "Programs for Digital Signal Processing", IEEE Press
%         John Wiley & Sons, 1979, Chap. 8.3.

% 20200203 Yuasa: modified from DECIMATE to apply matrix
% 20201214 Yuasa: ignore nan option
% 20220629 Yuasa: avoid error if data includes nan columns

narginchk(2,6);
error(nargoutchk(0,1,nargout,'struct'));

if nargin > 2
    nfilt  = convertStringsToChars(nfilt);
end

if nargin > 3
    option  = convertStringsToChars(option);
end

if nargin > 4
    dim     = convertStringsToChars(dim);
end

if nargin > 5
    nanflag = convertStringsToChars(nanflag);
end

% Validate required inputs 
validateinput(idata,r);

if fix(r) == 1
    odata = idata;
    return
end

fopt = 1;
empdim = true;
omitnan = false;
if nargin == 2
    nfilt = 8;
    dim  = 1;
else
    if ischar(nfilt)
        if nargin == 3 && (strcmpi(nfilt,'includenan')||strcmpi(nfilt,'omitnan'))
            nfilt = 8;
            dim  = 1;
            omitnan = strcmpi(nfilt,'omitnan');
        else
            if nfilt(1) == 'f' || nfilt(1) == 'F'
                fopt = 0;
            elseif nfilt(1) == 'i' || nfilt(1) == 'I'
                fopt = 1;
            else
                error(message('signal:decimate:InvalidEnum'))
            end
            if nargin > 5
                assert(strcmpi(nanflag,'includenan')||strcmpi(nanflag,'omitnan'),'Invalid option. Option must be ''omitnan'' or ''includenan''.');
                omitnan = strcmpi(nanflag,'omitnan');
            end
            if nargin > 4 && ~isempty(dim)
                empdim = false;
            else
                dim  = 1;
            end
            if nargin > 3 && ~isempty(option)
                nfilt = option;
            else
                nfilt = 8*fopt + 30*(1-fopt);
            end
        end
    else
        if nargin > 3
            if ischar(option)
                if nargin == 4 && (strcmpi(option,'includenan')||strcmpi(option,'omitnan'))
                    dim  = 1;
                    omitnan = strcmpi(option,'omitnan');
                    if isempty(nfilt)
                        nfilt = 8*fopt + 30*(1-fopt);
                    end
                else
                    if option(1) == 'f' || option(1) == 'F'
                        fopt = 0;
                    elseif option(1) == 'i' || option(1) == 'I'
                        fopt = 1;
                    else
                        error(message('signal:decimate:InvalidEnum'))
                    end
                    if nargin > 4 && ~isempty(dim)
                        empdim = false;
                    else
                        dim  = 1;
                    end
                    if isempty(nfilt)
                        nfilt = 8*fopt + 30*(1-fopt);
                    end
                end
            else
                if nargin > 4
                    if nargin > 5
                        assert(strcmpi(nanflag,'includenan')||strcmpi(nanflag,'omitnan'),'Invalid option. Option must be ''omitnan'' or ''includenan''.');
                        omitnan = strcmpi(nanflag,'omitnan');
                    end
                    if nargin == 5 && (strcmpi(dim,'includenan')||strcmpi(dim,'omitnan'))
                        omitnan = strcmpi(dim,'omitnan');
                    else
                        if dim(1) == 'f' || dim(1) == 'F'
                            fopt = 0;
                        elseif dim(1) == 'i' || dim(1) == 'I'
                            fopt = 1;
                        else
                            error(message('signal:decimate:InvalidEnum'))
                        end
                    end
                end
                if ~isempty(option)
                    dim    = option;
                    empdim = false;
                else
                    dim  = 1;
                end
                if isempty(nfilt)
                    nfilt = 8*fopt + 30*(1-fopt);
                end
            end
        else
            dim = 1;
            if isempty(nfilt)
                nfilt = 8*fopt + 30*(1-fopt);
            end
        end
    end
end
if fopt == 1 && nfilt > 13
    warning(message('signal:decimate:highorderIIRs'));
end

if isrow(idata) && empdim 	% row data
    dim = 2;
end
fulldim = max(ndims(idata),dim);
if omitnan
    idata = fillmissing(idata,'pchip',dim);
end
idata  = permute(idata,[dim setdiff(1:fulldim,dim)]);
datsiz = size(idata);
nd = datsiz(1);
nout = ceil(nd/r);
idata = reshape(idata,nd,[]);
nancol = all(isnan(idata),1);
idata(:,nancol) = 0;

if fopt == 0	% FIR filter
    b = fir1(nfilt,1/r);
    % prepend sequence with mirror image of first points so that transients
    % can be eliminated. then do same thing at right end of data sequence.
    nfilt = nfilt+1;
    itemp = 2*idata(1,:) - idata((nfilt+1):-1:2,:);
    [itemp,zi]=filter(b,1,itemp); %#ok
    [odata,zf] = filter(b,1,idata,zi);
    itemp = 2*idata(nd,:)-idata((nd-1):-1:(nd-2*nfilt),:);
    itemp = filter(b,1,itemp,zf);
    % finally, select only every r'th point from the interior of the lowpass
    % filtered sequence
    gd = grpdelay(b,1,8);
    list = round(gd(1)+1.25):r:nd;
    odata = odata(list,:);
    lod = length(list);
    nlen = nout - lod;
    nbeg = r - (nd - list(length(list)));
        odata = [odata; itemp(nbeg:r:nbeg+nlen*r-1,:)];
else	% IIR filter
    rip = .05;	% passband ripple in dB
    [b,a] = cheby1(nfilt, rip, .8/r);
    while all(b==0) || (abs(filtmag_db(b,a,.8/r)+rip)>1e-6)
        nfilt = nfilt - 1;
        if nfilt == 0
            break
        end
        [b,a] = cheby1(nfilt, rip, .8/r);
    end
    if nfilt == 0
        error(message('signal:decimate:InvalidRange'))
    end

    % be sure to filter in both directions to make sure the filtered data has zero phase
    % make a data vector properly pre- and ap- pended to filter forwards and back
    % so end effects can be obliterated.
    odata = filtfilt(b,a,idata);
    nbeg = r - (r*nout - nd);
    odata = odata(nbeg:r:nd,:);
end
odata(:,nancol) = nan;
odata = permute(reshape(odata,[size(odata,1),datsiz(2:end)]),[2:dim,1,(dim+1):fulldim]);

%--------------------------------------------------------------------------
function H = filtmag_db(b,a,f)
%FILTMAG_DB Find filter's magnitude response in decibels at given frequency.

nb = length(b);
na = length(a);
top = exp(-1i*(0:nb-1)*pi*f)*b(:);
bot = exp(-1i*(0:na-1)*pi*f)*a(:);

H = 20*log10(abs(top/bot));

%--------------------------------------------------------------------------
function validateinput(x,r)
% Validate 1st two input args: signal and decimation factor

if isempty(x) || issparse(x) || ~isa(x,'double')
    error(message('signal:decimate:invalidInput', 'X'));
end

if (abs(r-fix(r)) > eps) || (r <= 0)
    error(message('signal:decimate:invalidR', 'R'));
end

