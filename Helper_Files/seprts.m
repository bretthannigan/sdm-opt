function r = seprts(p)
% r = seprts(p)
% This program is for spectral factorization.
% The roots on the unit circle must have even degree.
% Roots with high multiplicity will cause problems,
% they should be handled by extracting them prior to
% using this program.
%
% See also: SFACT

SN = 0.0001; % Small Number (criterion for deciding if a
% root is on the unit circle).
rts = roots(p);

% The roots INSIDE the unit circle
irts = rts(abs(rts)<(1-SN));

% The roots ON the unit circle
orts = rts((abs(rts)>=(1-SN)) & (abs(rts)<=(1+SN)));
N = length(orts);
if rem(N,2) == 1
    disp('Sorry, but there is a problem (1) in seprts.m')
    r = [];
    return
end

% Sort roots on the unit circle by angle
[~,k] = sort(angle(orts));
orts = orts(k(1:2:end));

% Make final list of roots
r = [irts; orts];