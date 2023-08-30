function W = fweight(dcgain, wc, hfgain, order, varargin)
%FWEIGHT Create Weighting Transfer Function.
%   Arbitrary order weighting function with specified DC gain, crossover
%   frequency, and high frequency gain.
%
%   Inputs:
%       dcgain: magnitude as frequency approaches 0.
%       wc: crossover frequency (rad/s) where magnitude crosses 
%       N/sqrt(2) ~ N*3 dB, where N is the order.
%       hfgain: magnitude as frequency approaches Inf.
%       order: order of transfer function.
%       Ts: (optional) sample time (default=0).
%       dB: name-value parameter to indicate if dcgain and hfgain are
%       specified in decibels (default=false).
%       
%   Outputs:
%       W: LTI weight function in state-space form.
%   
%   See also: MAKEWEIGHT.M
%
%   $Author: BH $    $Date: 2017-10-27 $    $Revision: 1 $

%% Deal with variable function inputs.
p = inputParser;
addOptional(p, 'Ts', 0, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'dB', false, @(x) islogical(x) && isscalar(x));
parse(p, varargin{:});
isDecibels = p.Results.dB;
Ts = p.Results.Ts;

if isDecibels
    dcgain = db2mag(dcgain);
    hfgain = db2mag(hfgain);
end

z = zeros(1, order);
p = zeros(1, order);

%% Calculate poles and zeros.
p(:) = -wc;
z(:) = -wc*(dcgain/hfgain)^(1/order);
k = hfgain;

%% Generate LTI filter model.
W = zpk(z, p, k);
if Ts~=0
    W = c2d(W, Ts, 'matched');
end
W = ss(W); % state-space to maintain similar output to MAKEWEIGHT.M