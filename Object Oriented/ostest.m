function ostest(osr)

[b, a] = butter(6, 1/osr);

t = linspace(0, 1, 8192);
x1 = cos(2*pi*100*t);
x = filtfilt(b, a, x1);


xq1 = quantizer(x, 2^6);

t_os = linspace(0, 1, 8192*osr);
x_os = cos(2*pi*100*t_os);

xq21 = quantizer(x_os, 2^6);

xq22 = filtfilt(b, a, xq21);
xq2 = xq22(1:osr:end);
%xq2 = decimate(xq21, osr);

e1 = xq1 - x1;
e2 = xq2 - x1;

10*log10(var(x)/var(e1))
10*log10(var(x)/var(e2))

function y = quantizer(x, nlev)
    y = x*(nlev - 1);
    y = round(y)/(nlev - 1);
end
end