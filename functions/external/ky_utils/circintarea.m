function [S,rat1,rat2] = circintarea(x1,y1,r1,x2,y2,r2)

% S = CIRCINTAREA(x1,y1,r1,x2,y2,r2)
%   returns intersect area of two circles defined by center and radius in
%   x-y coordinates of <x1,y1,r1> and <x2,y2,r2>.
% 
% [S,ratio1,ratio2] = CIRCINTAREA(x1,y1,r1,x2,y2,r2)
%   also returns ratios of the intersect area by each area of the circles.

% 20210803 Yuasa

narginchk(6,6);

%-- Sub-parameters
%%-- d^2  = (x1-x2)^2 + (y1-y2)^2
%%-- r1^2 = x^2 + y^2
%%-- r2^2 = (d-x)^2 + y^2
d = sqrt((x1-x2).^2+(y1-y2).^2);
x = (r1.^2 - r2.^2 + d.^2)./(2.*d);
y = sqrt(r1.^2 - x.^2);

%-- Outputs
%%-- S_intetsect = S_sector1 + S_sector2 - S_rhombus
%%-- S_sector1   = r1^2 * pi * theta1 / 2pi
%%-- S_sector2   = r2^2 * pi * theta2 / 2pi
%%-- S_rhombus   = 2 * d * y / 2
%%-- theta1 = acos(x/r1)
%%-- theta2 = acos((d-x)/r2)
S = real(r1.^2 .* acos(x./r1) + r2.^2 .* acos((d-x)./r2) - d.*y);
rat1 = S ./ (r1.^2 * pi);
rat2 = S ./ (r2.^2 * pi);