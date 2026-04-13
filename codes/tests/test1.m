%% Experiment 1: Scaling law for linear CLF-QP  control
% Figure (a) and (b)

clear; close all; clc;

%% ----------------------- User settings -----------------------
dt    = 1e-2;      % [s]
sim_t = 5.0;       % [s]

% Parameters
params.m = 1;
params.l = 1;
params.g = 9.81;
params.b = 0.01;
params.I = params.m*params.l^2/3;   % uniform rod inertia about pivot

params.u_max = 20;
params.u_min = -params.u_max;

% CLF design
params.Kp = 6;
params.Kd = 5;

% experiment sweep
sigma_vals = [0.1, 0.3, 0.5, 1, 1.5, 2, 3];
psi0_vals  = [75, 60, 45, 30, 15] * pi/180;
omega0     = 0.0;

% mini-norm hard CLF-QP
with_slack_linear = false;

%% ----------------------- Plot styling -----------------------
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesFontSize',16);
set(groot,'defaultAxesLineWidth',1.2);
set(groot,'defaultAxesFontName','Times New Roman');
set(groot,'defaultTextFontName','Times New Roman');

%% ----------------------- Build CLF -----------------------
[P, ~, ~] = build_quadratic_clf(params);

%% ----------------------- Sweep -----------------------
nPsi   = numel(psi0_vals);
nSigma = numel(sigma_vals);

Uinf_mat       = zeros(nPsi, nSigma);
MeanSolve      = nan(nPsi, nSigma);
Infeasible_mat = false(nPsi, nSigma);

for ip = 1:nPsi
    psi0 = psi0_vals(ip);
    x0   = [psi0; omega0];

    for is = 1:nSigma
        params.sigma = sigma_vals(is);

        traj = simulate_zoh_linear(x0, dt, sim_t, params.u_max, params, P, with_slack_linear);

        Uinf_mat(ip,is) = max(abs(traj.u));

        st = traj.solve_t(~isnan(traj.solve_t));
        if ~isempty(st)
            MeanSolve(ip,is) = mean(st);
        end

        if any(isnan(traj.slack))
            Infeasible_mat(ip,is) = true;
        end
    end
end

%% ----------------------- Common colors/styles -----------------------
% red, black, blue, purple, green
cols_phi = [
    1.0000, 0.0000, 0.0000;   % red
    0.0000, 0.0000, 0.0000;   % black
    0.0000, 0.4470, 0.7410;   % blue
    0.4940, 0.1840, 0.5560;   % purple
    0.4660, 0.6740, 0.1880    % green
];

mks   = {'o-','s-','d-','^-','v-'};
ls_ic = {'-','--',':','-.','-'};

%% ----------------------- Figure 1(a): ||u||_inf vs sigma -----------------------
figure(1); clf;
set(gcf,'Color','w','Name','Figure 1');
hold on; grid on; box on;

for ip = 1:nPsi
    plot(sigma_vals, Uinf_mat(ip,:), mks{ip}, ...
        'Color', cols_phi(ip,:), ...
        'LineWidth', 2.2, ...
        'MarkerSize', 8, ...
        'DisplayName', sprintf('$\\phi_0=%d^\\circ$', round(psi0_vals(ip)*180/pi)));
end

for ip = 1:nPsi
    infeas_idx = find(Infeasible_mat(ip,:));
    if ~isempty(infeas_idx)
        plot(sigma_vals(infeas_idx), Uinf_mat(ip,infeas_idx), 'x', ...
            'Color', cols_phi(ip,:), ...
            'LineWidth', 2.0, ...
            'MarkerSize', 10, ...
            'HandleVisibility','off');
    end
end

xlabel('$\sigma$','FontSize',22);
ylabel('$\|u\|_{\infty}$','FontSize',22);
ylim([0 16]);
legend('Location','northwest');

%% ----------------------- Figure 1(b): u versus V for all initial angles and sigmas -----------------------
figure(2); clf;
set(gcf,'Color','w','Name','Figure 4');
hold on; grid on; box on;

% color -> sigma
cols_sigma = [
    0.8500, 0.0000, 0.0000;   % red
    0.0000, 0.0000, 0.0000;   % black
    0.0000, 0.4470, 0.7410;   % blue
    0.4940, 0.1840, 0.5560;   % purple
    0.4660, 0.6740, 0.1880;   % green
    0.6350, 0.0780, 0.1840;   % dark red
    0.3010, 0.7450, 0.9330    % cyan
];

% line style -> initial value
ls_ic = {'-','--',':','-.','-'};

% maximum Lyapunov level over all initial conditions
Vmax_all = 0;
for ip = 1:nPsi
    x0_tmp = [psi0_vals(ip); omega0];
    V0_tmp = x0_tmp.' * P * x0_tmp;
    Vmax_all = max(Vmax_all, V0_tmp);
end

% store one anchor point for each initial value, used for phi0 label
label_x = zeros(nPsi,1);
label_y = zeros(nPsi,1);

for ip = 1:nPsi
    psi0 = psi0_vals(ip);
    x0   = [psi0; omega0];

    for is = 1:nSigma
        params.sigma = sigma_vals(is);
        traj = simulate_zoh_linear(x0, dt, sim_t, params.u_max, params, P, with_slack_linear);

        plot(traj.V, traj.u, ...
            'LineStyle', ls_ic{ip}, ...
            'Color', cols_sigma(is,:), ...
            'LineWidth', 1.6, ...
            'HandleVisibility','off');

        % mark start point
        plot(traj.V(1), traj.u(1), 'o', ...
            'Color', cols_sigma(is,:), ...
            'MarkerFaceColor', cols_sigma(is,:), ...
            'MarkerSize', 4.5, ...
            'HandleVisibility','off');

        % use one representative curve to place the phi0 label
        % here choose the middle sigma for label anchor
        if is == ceil(nSigma/2)
            idx_txt = max(3, round(0.18*numel(traj.V)));
            label_x(ip) = traj.V(idx_txt);
            label_y(ip) = traj.u(idx_txt);
        end
    end
end

% phi0 labels: black, corresponding to different initial values
for ip = 1:nPsi
    text(label_x(ip), label_y(ip), ...
        sprintf('$\\phi_0=%d^\\circ$', round(psi0_vals(ip)*180/pi)), ...
        'Color', 'k', ...
        'FontSize', 15, ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','bottom');
end

% legend for sigma only
h_sigma = gobjects(1,nSigma);
for is = 1:nSigma
    h_sigma(is) = plot(nan, nan, '-', ...
        'Color', cols_sigma(is,:), ...
        'LineWidth', 2.0, ...
        'DisplayName', sprintf('$\\sigma=%.1f$', sigma_vals(is)));
end

xlabel('$V(x)$','FontSize',22);
ylabel('$u$','FontSize',22);
xlim([0 Vmax_all]);
legend(h_sigma, 'Location','best');
