function ordstr = int2ordinal(number)

% INT2ORDINAL convert positive integer to ordinal character vector.
% 
% Usage:
%   ordstr = INT2ORDINAL(number)
% 
% Example: 
% The following example returns the character vector '4th'.
%   str = int2ordinal(4)
% 
% The following example returns the character vector 'fourth'.
%   str = iptnum2ordinal(4)
% 
% See also: IPTNUM2ORDINAL

number = round(abs(number(1)));

if number == 0
    ordstr = '0th';
else
    ordstr = iptnum2ordinal(number);
    ordstr = [int2str(number) ordstr((end-1):end)];
end