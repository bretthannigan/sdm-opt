function [h,r] = sfact(p)
% [h,r] = sfact(p)
% spectral factorization of a polynomial p.
% h: new polynomial
% r: roots of new polynomial
%
% % example:
%   g = rand(1,10);
%   p = conv(g,g(10:-1:1));
%   h = sfact(p);
%   p - conv(h,h(10:-1:1)) % should be 0
% See also: SEPRTS, LEJA
% leja.m is by Markus Lang, and is available from the
% Rice University DSP webpage: http://www.dsp.rice.edu/

if length(p) == 1
    h = p;
    r = [];
    return
end

% Get the appropriate roots.
r = seprts(p);

% Form the polynomial from the roots
r = leja(r);
h = poly(r);
if isreal(p)
    h = real(h);
end

% normalize
h = h*sqrt(max(p)/sum(abs(h).^2));