function H_loop = sigma_delta_ltf(H, a, b)
%SIGMA_DELTA_LTF Sigma Delta Modulator Noise Transfer Function
%   The noise transfer function of an even-order sigma-delta modulator with
%   multiple feedback chain of loop filters topology. 
%
%   INPUTS:
%       H: filter transfer function object. If vector, transfer functions 
%       are different. If scalar, all the chain transfer functions are 
%       identical.
%       a: vector of m feedback coefficients, where 2m is the order of the
%       sigma delta modulator.
%       b: vector of n integrator weights, where n is the order of the
%       sigma delta modulator.
%
%   OUTPUT:
%       H_loop: total loop transfer function object.
%
%   See also: 
%       [1] T. Ritoniemi, T. Krema, H. Tenhunen, "Design of stable High 
%           Order 1-bit Sigma-Delta Modulators", IEEE Proc. IS-CAS'90, 
%           pp. 3267-3270, 1990.

    B = 0;  A = 1; H_loop = 1;
    n = length(b);
    m = length(a);
    if ~isrow(a)
        a = a.';
    end
    if ~isrow(b)
        b = b.';
    end
    if mod(n, 2)~=0
        error('sigma_delta_ltf:oddOrder', 'The order n=length(b) of the sigma delta modulator must be even.')
    end
    if 2*m~=n
        error('sigma_delta_ltf:invalidm', ['The length of the a vector must be equal to ' num2str(n/2) '.'])
    end
    if isscalar(H)
        H = H*ones(n, 1);
    else
        if length(H)~=n
            error('sigma_delta_ltf:invalidH', ['H must be a scalar or a vector with length ' num2str(n) '.'])
        end
    end

    for ii=1:m
        B_jj = 1;
        for jj=ii+1:m
            B_jj = B_jj*(1 + a(jj)*minreal(H(2*jj-1)*H(2*ii)));
        end
        B_kk = 1;
        for kk=1:ii-1
            B_kk = B_kk*minreal(H(2*kk-1)*H(2*kk));
        end
        B = B + minreal(b(2*ii-1)*H(2*ii-1) + b(2*ii)*minreal(H(2*ii-1)*H(2*ii)))*B_jj*B_kk;
        A = A*(1 + a(ii)*minreal(H(2*ii-1)*H(2*ii)));
        H_loop = H_loop*(B/A);
    end
end

