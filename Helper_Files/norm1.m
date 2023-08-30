function n = norm1(sys)
%NORM1 Approximate l1 Norm of LTI System
%   Calculates the (approximate) l1 norm based on the impulse response of
%   the LTI object sys. Output is approximate because the impulse response
%   is terminated at a point decided by IMPULSE.M.
%
%   Input:
%       sys: linear time invariant system object.
%   Output:
%       n: p x q numeric matrix of the l1 norm from input p to output q of
%       sys.
%
%   See also: NORM, IMPULSE
%
%   $Author: BH $    $Date: 2018-01-25 $    $Revision: 1 $
    
    if ~isa(sys, 'lti')
        error('norm1:notLti', 'sys must be an LTI object.')
    end     
    isInfinite = ~arrayfun(@isstable, sys);
    if sys.Ts==0 % Continuous time - integrate the impulse response.
        [y, t] = impulse(sys);
        n = trapz(t, y);
    else % Discrete time - sum the impulse response.
        n = sum(abs(impulse(sys)*sys.Ts));
    end
    n = reshape(squeeze(n), size(sys));
    n(isInfinite) = Inf;
end