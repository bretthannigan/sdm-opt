% Equations to reproduce Figure 6.4 from [Risbo1994]

H2_gauss = @(mu_y) (1 - mu_y.^2)./(1 - mu_y.^2 - (2/pi)*exp(-2*(erfinv(mu_y).^2)));
H2_uniform = @(mu_y) (1 - mu_y.^2)./(1 - mu_y.^2 - (3/4)*(1 - mu_y.^2).^2);
H2_tri = @(mu_y) (1 - mu_y.^2)./(1 - mu_y.^2 - (2/3)*(3*(1 - mu_y) - 2*(1 - mu_y).^(3/2)).^2);
mu_y = 0:0.001:1;
figure
plot(mu_y, H2_gauss(mu_y))
hold on
plot(mu_y, H2_uniform(mu_y))
plot(mu_y, H2_tri(mu_y))

H2inv_gauss = @(H2) fminsearch(@(x) abs(H2_gauss(x)-H2), 0);
H2inv_uniform = @(H2) fminsearch(@(x) abs(H2_uniform(x)-H2), 0);
H2inv_tri = @(H2) fminsearch(@(x) abs(H2_tri(x)-H2), 0);