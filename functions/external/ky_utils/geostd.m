function varargout = geostd(x,varargin)

%GEOSTD Geometric standard deviation.
%   Y = GEOSTD(X) returns the geometric standard deviation of the values in
%   X.  For matrices, Y is a row vector containing the geometric standard
%   deviation of each column. For N-D arrays, GEOSTD operates along the first
%   non-singleton dimension of X.
%
%   GEOSTD normalizes Y by N-1 if N>1, where N is the sample size.  This is
%   the sqrt of an unbiased estimator of the variance of the population
%   from which X is drawn, as long as X consists of independent,
%   identically distributed samples. For N=1, Y is normalized by N.
%
%   Y = GEOSTD(X,1) normalizes by N and produces the square root of the second
%   moment of the sample about its mean.  GEOSTD(X,0) is the same as GEOSTD(X).
%
%   Y = GEOSTD(X,W) computes the standard deviation using the weight vector W.
%   W typically contains either counts or inverse variances.  The length of
%   W must equal the length of the dimension over which GEOSTD operates, and
%   its elements must be nonnegative.  If X(I) is assumed to have standard
%   deviation proportional to 1/SQRT(W(I)), then Y * SQRT(MEAN(W)/W(I)) is
%   an estimate of the standard deviation of X(I).  In other words, Y *
%   SQRT(MEAN(W)) is an estimate of standard deviation for an observation
%   given weight 1.
%
%   Y = GEOSTD(X,0,'all') or Y = GEOSTD(X,1,'all') returns the geometric
%   standard deviation of all elements of X. A weight of 0 normalizes by
%   N-1 and a weight of 1 normalizes by N.
%
%   Y = GEOSTD(X,W,DIM) takes the standard deviation along the dimension DIM
%   of X.
%
%   Y = GEOSTD(X,0,VECDIM) or Y = GEOSTD(X,1,VECDIM) operates on the dimensions 
%   specified in the vector VECDIM. A weight of 0 normalizes by N-1 and a 
%   weight of 1 normalizes by N. For example, GEOSTD(X,0,[1 2]) operates on
%   the elements contained in the first and second dimensions of X.
%
%   The standard deviation is the square root of the variance (VAR).
%
%   GEOSTD(...,NANFLAG) specifies how NaN (Not-A-Number) values are treated.
%   The default is 'includenan':
%
%   'includenan' - the standard deviation of a vector containing NaN 
%                  values is also NaN.
%   'omitnan'    - elements of X or W containing NaN values are ignored.
%                  If all elements are NaN, the result is NaN.
% 
%   [Yn Yp] = GEOSTD(X,...) returns the size of the both side of geometric
%   standard deviation separately.
% 
%   See also STD, GEOMEAN.

%   20211217 Yuasa

if any(x(:) < 0) || ~isreal(x)
    error(message('stats:geomean:BadData'))
end
nargoutchk(0,2);

% Take the n-th root of the product of elements of X, along dimension DIM.
y = exp(std(log(x),varargin{:}));

if nargout <= 1
    varargout{1} = y;
else
    m = exp(mean(log(x),varargin{2:end}));
    yn = exp(log(m) - log(y));
    yp = exp(log(m) + log(y));
    varargout{1} = m - yn;
    varargout{2} = yp - m;
end

