function yy = modmod(xx,mm)

% MODMOD   modulo operation
% 
%   y = MODMOD(x,m)
% 
%   This function outputs 'm' instead of '0',
%   when x is a multiple of m.
% 

% 20160203 Yuasa
% 20160525 Bug fix
% 20200331 update help document

yy = mod(xx-1,mm)+1;