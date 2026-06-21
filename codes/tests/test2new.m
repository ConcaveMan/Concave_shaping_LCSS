%% Experiment 2: Linear vs concave-shaped CLF-QP
% Figure (c) and (d)
% Fixed umax = 10, sigma = 3
% Concave-shaped CLF-QP with kmin = 0.1, kmax = 2, r = [1,0.9,0.8,0.7,0.6]
% Plot normalized V(t) and u-V, and collect energy / nominal rate / convergence time

clear; close all; clc;

%% ----------------------- User settings -----------------------
dt    = 1e-2;      % [s]
sim_t = 5.0;       % [s]

% Parameters
params.m = 1;
params.l = 1;
params.g = 9.81;
params.b = 0.01;
params.I = params.m*params.l^2/3;

params.u_max = 10;
params.u_min = -params.u_max;

params.Kp = 6;
params.Kd = 5;

params.sigma = 3;

with_slack_linear  = false;
with_slack_concave = false;

% initial condition
x0 = [pi/4; 0.0];

% concave factor settings
kmin = 0.1;
kmax = 2.2;
p    = 1.0;
r_vals = [1.0, 0.9, 0.8, 0.7, 0.6];

%% ----------------------- Plot styling -----------------------
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesFontSize',16);
set(groot,'defaultAxesLineWidth',1.2);
set(groot,'defaultAxesFontName','Times New Roman');
set(groot,'defaultTextFontName','Times New Roman');

%% ----------------------- Build CLF -----------------------
[P, Vfun, ~] = build_quadratic_clf(params);
c = Vfun(x0);

eps_levels = [1e-1, 1e-2, 1e-3] * c;
eps_ratio_vals = [1e-1; 1e-2; 1e-3];

%% ----------------------- Simulate linear -----------------------
N  = ceil(sim_t/dt) + 1;
nx = numel(x0);

X = zeros(N,nx);
t = zeros(N,1);
u = zeros(N-1,1);
V = zeros(N-1,1);
s = zeros(N-1,1);
rho = nan(N-1,1);
solve_t = nan(N-1,1);
tu = zeros(N-1,1);

X(1,:) = x0.';
t(1)   = 0;

xk = x0;
tk = 0;

for k = 1:N-1
    [uk, sk, Vk, rhok, tk_solve] = clf_qp_linear(xk, params, P, with_slack_linear);

    if isempty(uk), uk = 0; else, uk = uk(1); end
    if isempty(sk), sk = NaN; else, sk = sk(1); end
    if isempty(Vk), Vk = xk.' * P * xk; else, Vk = Vk(1); end
    if isempty(rhok), rhok = nan; else, rhok = rhok(1); end
    if isempty(tk_solve), tk_solve = nan; else, tk_solve = tk_solve(1); end

    uk = min(max(uk, params.u_min), params.u_max);

    tu(k)      = tk;
    u(k)       = uk;
    V(k)       = Vk;
    s(k)       = sk;
    rho(k)     = rhok;
    solve_t(k) = tk_solve;

    ode = @(tt,xx) pendulum_dyn(xx, uk, params);
    [~, xs] = ode45(ode, [0 dt], xk);

    xk = xs(end,:).';
    tk = tk + dt;

    X(k+1,:) = xk.';
    t(k+1)   = tk;
end

traj_lin.t  = t;
traj_lin.X  = X;
traj_lin.tu = tu;
traj_lin.u  = u;
traj_lin.V  = V;
traj_lin.slack = s;
traj_lin.rho = rho;
traj_lin.solve_t = solve_t;

%% ----------------------- Simulate concave family -----------------------
nr = numel(r_vals);
R = struct('r',[],'traj',[],'label','');

for ir = 1:nr
    r_i = r_vals(ir);

    N  = ceil(sim_t/dt) + 1;
    nx = numel(x0);

    X = zeros(N,nx);
    t = zeros(N,1);
    u = zeros(N-1,1);
    V = zeros(N-1,1);
    s = zeros(N-1,1);
    rho = nan(N-1,1);
    solve_t = nan(N-1,1);
    tu = zeros(N-1,1);

    X(1,:) = x0.';
    t(1)   = 0;

    xk = x0;
    tk = 0;

    for k = 1:N-1
        [uk, sk, Vk, rhok, tk_solve] = clf_qp_concave(xk, params, P, c, r_i, kmin, kmax, p, with_slack_concave);

        if isempty(uk), uk = 0; else, uk = uk(1); end
        if isempty(sk), sk = NaN; else, sk = sk(1); end
        if isempty(Vk), Vk = xk.' * P * xk; else, Vk = Vk(1); end
        if isempty(rhok), rhok = nan; else, rhok = rhok(1); end
        if isempty(tk_solve), tk_solve = nan; else, tk_solve = tk_solve(1); end

        uk = min(max(uk, params.u_min), params.u_max);

        tu(k)      = tk;
        u(k)       = uk;
        V(k)       = Vk;
        s(k)       = sk;
        rho(k)     = rhok;
        solve_t(k) = tk_solve;

        ode = @(tt,xx) pendulum_dyn(xx, uk, params);
        [~, xs] = ode45(ode, [0 dt], xk);

        xk = xs(end,:).';
        tk = tk + dt;

        X(k+1,:) = xk.';
        t(k+1)   = tk;
    end

    traj_i.t  = t;
    traj_i.X  = X;
    traj_i.tu = tu;
    traj_i.u  = u;
    traj_i.V  = V;
    traj_i.slack = s;
    traj_i.rho = rho;
    traj_i.solve_t = solve_t;

    R(ir).r = r_i;
    R(ir).traj = traj_i;
    R(ir).label = sprintf('$r=%.1f$', r_i);
end

%% ----------------------- Common colors -----------------------
cols = [
    0.0000, 0.0000, 0.0000;   % linear black
    1.0000, 0.0000, 0.0000;   % red
    0.0000, 0.4470, 0.7410;   % blue
    0.4940, 0.1840, 0.5560;   % purple
    0.4660, 0.6740, 0.1880;   % green
    0.8500, 0.3250, 0.0980    % orange
];

%% ----------------------- Plot (c): normalized V(t) -----------------------
figure('Color','w','Name','Experiment 2: normalized V(t)');
hold on; grid on; box on;
set(gca,'YScale','log');

V0 = c;

plot(traj_lin.tu, max(traj_lin.V./V0, 1e-12), '-', ...
    'Color', cols(1,:), 'LineWidth', 2.4, 'DisplayName', 'linear');

for ir = 1:nr
    plot(R(ir).traj.tu, max(R(ir).traj.V./V0, 1e-12), '-', ...
        'Color', cols(ir+1,:), ...
        'LineWidth', 2.2, ...
        'DisplayName', R(ir).label);
end

% epsilon lines (not shown in legend)
yline(1e-1, 'k--', 'LineWidth', 1.2, 'HandleVisibility','off');
yline(1e-2, 'k-.', 'LineWidth', 1.2, 'HandleVisibility','off');
yline(1e-3, 'k:',  'LineWidth', 1.4, 'HandleVisibility','off');

xlabel('$t$','FontSize',22);
ylabel('$V(x(t))$','FontSize',22);

ylim([1e-4 1]);
yticks([1e-4 1e-3 1e-2 1e-1 1e0]);
yticklabels({'$10^{-4}V_0$','$10^{-3}V_0$','$10^{-2}V_0$','$10^{-1}V_0$','$V_0$'});

% put epsilon labels on the lines
x_txt = 0.46 * max(traj_lin.tu);

text(x_txt, 1e-1, '$c=V_0, ~\epsilon=0.1c$', ...
    'FontSize', 15, ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','bottom', ...
    'Color', 'k');

text(x_txt, 1e-2, '$c=V_0, ~\epsilon=0.01c$', ...
    'FontSize', 15, ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','bottom', ...
    'Color', 'k');

text(x_txt, 1e-3, '$c=V_0, ~\epsilon=0.001c$', ...
    'FontSize', 15, ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','bottom', ...
    'Color', 'k');

lgd = legend('Location','northeast');
lgd.Box = 'on';
%% ----------------------- Plot (d): u versus V -----------------------
figure('Color','w','Name','Experiment 2: u versus V');
hold on; grid off; box on;

% linear peak input
upeak_lin = max(abs(traj_lin.u));

plot(traj_lin.V, traj_lin.u, '-', ...
    'Color', cols(1,:), 'LineWidth', 2.4, ...
    'DisplayName', sprintf('linear, $\\|u\\|_{\\infty}=%.2f$', upeak_lin));

for ir = 1:nr
    upeak_i = max(abs(R(ir).traj.u));
    peak_reduction = upeak_lin - upeak_i;
    reduction_pct = 100 * peak_reduction / max(upeak_lin,1e-12);

    plot(R(ir).traj.V, R(ir).traj.u, '-', ...
        'Color', cols(ir+1,:), ...
        'LineWidth', 2.2, ...
        'DisplayName', sprintf('$r=%.1f$, $\\|u\\|_{\\infty}=%.2f$, reduced peak $(%.1f\\%%)$', ...
        R(ir).r, upeak_i, reduction_pct));
end

xlabel('$V(x)$','FontSize',22);
ylabel('$u$','FontSize',22);
xlim([0 c]);
legend('Location','best');
legend boxoff

%% ----------------------- Plot 3: u versus t -----------------------
figure('Color','w','Name','Experiment 2: u versus t');
hold on; grid on; box on;

plot(traj_lin.tu, traj_lin.u, '-', ...
    'Color', cols(1,:), 'LineWidth', 2.4, 'DisplayName', 'linear');

for ir = 1:nr
    plot(R(ir).traj.tu, R(ir).traj.u, '-', ...
        'Color', cols(ir+1,:), ...
        'LineWidth', 2.2, ...
        'DisplayName', R(ir).label);
end

xlabel('$t$','FontSize',22);
ylabel('$u(t)$','FontSize',22);
legend('Location','best');

%% ----------------------- Plot 4: sigma_alpha versus epsilon/c -----------------------
figure('Color','w','Name','Experiment 2: windowed nominal rate');
hold on; grid on; box on;

eps_grid = logspace(-4, -0.02, 400);   % epsilon/c
sigma_base = params.sigma;

% baseline
yline(sigma_base, 'k', 'LineWidth', 3, ...
    'DisplayName', '$\sigma=3$');

for ir = 1:nr
    r_i = R(ir).r;
    sig_curve = zeros(size(eps_grid));

    for jj = 1:numel(eps_grid)
        eps_j = eps_grid(jj)*c;
        sig_curve(jj) = rational_window_rate(eps_j,c,params.sigma,kmin,kmax,r_i);
    end

    plot(eps_grid, sig_curve, '-', ...
        'Color', cols(ir+1,:), ...
        'LineWidth', 2.2, ...
        'DisplayName', sprintf('$r=%.1f$', r_i));
end

% mark the three windows used in simulations
xline(1e-1, 'k--', 'LineWidth', 1.0, 'HandleVisibility','off');
xline(1e-2, 'k-.', 'LineWidth', 1.0, 'HandleVisibility','off');
xline(1e-3, 'k:',  'LineWidth', 1.0, 'HandleVisibility','off');

xlabel('$\epsilon/c$','FontSize',22);
ylabel('$\sigma_{\alpha}(\epsilon,c)$','FontSize',22);

set(gca,'XScale','log');
xlim([1e-4 1]);
legend('Location','best');
legend boxon

%% ----------------------- Results table -----------------------
Controller = {};
EpsRatio = [];
Tconv = [];
Energy = [];
NominalRate = [];

% linear
for ie = 1:numel(eps_levels)
    epsV = eps_levels(ie);
    idx = find(traj_lin.V <= epsV, 1, 'first');
    if isempty(idx), idx = numel(traj_lin.V); end

    Tc = traj_lin.tu(idx);
    E  = sum(traj_lin.u(1:idx).^2) * dt;

    V0i = traj_lin.V(1);
    Ve = traj_lin.V(idx);
    sigma_nom = (1/max(Tc,1e-12)) * log(max(V0i,1e-12)/max(Ve,1e-12));

    Controller{end+1,1} = 'linear';
    EpsRatio(end+1,1)   = eps_ratio_vals(ie);
    Tconv(end+1,1)      = Tc;
    Energy(end+1,1)     = E;
    NominalRate(end+1,1)= sigma_nom;
end

% concave family
for ir = 1:numel(R)
    for ie = 1:numel(eps_levels)
        epsV = eps_levels(ie);
        idx = find(R(ir).traj.V <= epsV, 1, 'first');
        if isempty(idx), idx = numel(R(ir).traj.V); end

        Tc = R(ir).traj.tu(idx);
        E  = sum(R(ir).traj.u(1:idx).^2) * dt;

        V0i = R(ir).traj.V(1);
        Ve = R(ir).traj.V(idx);
        sigma_nom = (1/max(Tc,1e-12)) * log(max(V0i,1e-12)/max(Ve,1e-12));

        Controller{end+1,1} = sprintf('r=%.1f', R(ir).r);
        EpsRatio(end+1,1)   = eps_ratio_vals(ie);
        Tconv(end+1,1)      = Tc;
        Energy(end+1,1)     = E;
        NominalRate(end+1,1)= sigma_nom;
    end
end

Results = table(Controller, EpsRatio, Tconv, Energy, NominalRate);

controller_order = [{'linear'}, arrayfun(@(rr) sprintf('r=%.1f', rr), r_vals, 'UniformOutput', false)];
Results.Controller = categorical(Results.Controller, controller_order, 'Ordinal', true);
Results = sortrows(Results, {'Controller','EpsRatio'}, {'ascend','descend'});

disp('--- Experiment 2: energy / nominal rate / convergence time ---');
disp(Results);

function sig = rational_window_rate(eps,c,sigmaL,kmin,kmax,r)
    ell = (r-kmin)*c/(kmax-r);
    T = (1/sigmaL)*( ...
        (1/kmax)*log(c/eps) + ...
        (kmax-kmin)/(kmax*kmin)* ...
        log((kmin*c+kmax*ell)/(kmin*eps+kmax*ell)) );
    sig = log(c/eps)/T;
end