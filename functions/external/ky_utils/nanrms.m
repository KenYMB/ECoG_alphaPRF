function y = nanrms(x,dim)
%NANRMS Root mean squared value, ignoring NaNs.
%   For vectors, NANRMS(X) is the root mean squared value of the non-NaN
%   elements in X. For matrices, RMS(X) is a row vector containing the RMS
%   value of the non-NaN elements from each column. For N-D arrays,
%   NANRMS(X) operates along the first non-singleton dimension.
%
%   Y = NANRMS(X,DIM) operates along the dimension DIM.
%
%   When X is complex, the NANRMS is computed using the magnitude
%   NANRMS(ABS(X)). 
%
%   See also RMS, NANMEAN, NANMEDIAN, NANSTD, NANVAR, NANMIN, NANMAX, NANSUM.

%   20181210 Yuasa

if nargin==1
  y = sqrt(nanmean(x .* conj(x)));
else
  y = sqrt(nanmean(x .* conj(x), dim));
end
