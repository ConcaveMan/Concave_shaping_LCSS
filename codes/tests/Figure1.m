clear; clc; close all;

sigmaL = 1;
c = 10;
eps_list = [1, 0.1, 0.01];

kmax = 5;
kmin = 0.1;

gamma = kmin*(kmax-1)/(kmax-kmin);

v = linspace(0,c,1000);
alpha_l = sigmaL*v;

figure; hold on; grid on; box on;
plot(v, alpha_l, 'k--', 'LineWidth', 2, 'DisplayName', '$\sigma_L v$');

for eps = eps_list

    % Conservative sufficient threshold:
    % eps/c <= (beta/(1+beta))^(1/gamma)
    a = (eps/c)^gamma;
    beta_req = a/(1-a);

    % beta(r)= kmax*(r-kmin)/(kmin*(kmax-r))
    r_bound = (beta_req*kmin*kmax + kmax*kmin) / ...
              (kmax + beta_req*kmin);

    % Exact break-even r from sigma_alpha(eps,c)=sigmaL
    rate_fun = @(r) window_rate(eps,c,sigmaL,kmin,kmax,r) - sigmaL;
    r_exact = fzero(rate_fun,[kmin+1e-8,1]);

    % Use conservative r for plotting
    r = r_bound;
    ell = (r-kmin)*c/(kmax-r);
    s = (kmin*v + kmax*ell)./(v + ell);
    alpha_cave = sigmaL*s.*v;

    sig_val = window_rate(eps,c,sigmaL,kmin,kmax,r);

    plot(v, alpha_cave, 'LineWidth', 2, ...
        'DisplayName', sprintf('$\\epsilon/c=%.0e$, $r=%.2f$, $\\sigma_\\alpha(\\epsilon,c)=%.2f$', ...
        eps/c, r, sig_val));

    fprintf('eps = %.4g, eps/c = %.4g\n', eps, eps/c);
    fprintf('  conservative r_min = %.6f, sigma_alpha = %.6f\n', r_bound, sig_val);
    fprintf('  exact break-even r = %.6f\n\n', r_exact);
end

xlabel('$v$','Interpreter','latex');
ylabel('$\alpha(v)$','Interpreter','latex');
legend('Interpreter','latex','Location','northwest');
title('Endpoint-relaxed concave comparison laws','Interpreter','latex');


function sig = window_rate(eps,c,sigmaL,kmin,kmax,r)
    ell = (r-kmin)*c/(kmax-r);
    T = (1/sigmaL)*( ...
        (1/kmax)*log(c/eps) + ...
        (kmax-kmin)/(kmax*kmin)* ...
        log((kmin*c+kmax*ell)/(kmin*eps+kmax*ell)) );
    sig = log(c/eps)/T;
end
