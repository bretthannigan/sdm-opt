function [y, f] = fft_handler(varargin)

    %% Deal with variable function inputs.
    p = inputParser;
    addRequired(p, 'x', @(x) isnumeric(x) && isvector(x));
    parse(p, varargin{1});
    x = p.Results.x; % Parse to obtain length(x).
    % If the sampling frequency is specified:
    addOptional(p, 'f_s', 2*pi, @(x) isnumeric(x) && isscalar(x));
    % If a windowing function is specified:
    addParameter(p, 'window', rectwin(length(x)), @(x) isnumeric(x) && isvector(x));
    parse(p, varargin{:});
    f_s = p.Results.f_s;
    w = p.Results.window;
    
    %% Calculate signal lengths.
    x = w.*x;
    nfft = 2^nextpow2(length(x));
    y = fft(x, nfft);
    
    ny = ceil((nfft+1)/2);
    f = linspace(0, f_s/2, ny);
    %y = abs(y(1:ny));
    y = 2*abs(y(1:ny))/length(y);
    
    %% Plot if necessary.
    if nargout==0
        %figure;
        plot(2*pi*f, 20*log10(y));
        title('Frequency Spectrum')
        if ismember('f_s', p.UsingDefaults) % Normalized frequency.
            xlabel('Normalized frequency (rad/sample)')
        else
            xlabel('Frequency (rad/s)')
        end
        ylabel('Magnitude (dB)')
        grid on;
    end
end