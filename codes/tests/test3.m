%% Experiment 3: linear vs optimal-decay vs concave-shaped CLF-QP
% x0 = [pi/4; 0], theta = 10
% Compare:
%   1) linear CLF-QP with sigma = 3
%   2) optimal-decay CLF-QP with q_delta = 1e5, q_rho = 100
%   3) optimal-decay CLF-QP with q_delta = 1e5, q_rho = 200
%   4) optimal-decay CLF-QP with q_delta = 1e5, q_rho = 300
%   5) optimal-decay CLF-QP with q_delta = 1e5, q_rho = 500
%   6) concave-shaped CLF-QP with baseline sigma = 3 and r = 1
%
% Plots:
%   Figure (e): normalized V(t)
%   Figure (f): u versus V, legend shows total energy integral
%
% Metrics:
%   epsilon/c = 1e-1, 1e-2, 1e-3
%   convergence time, energy, nominal rate, peak input, mean solve time, integrated slack

clear; close all; clc;

%% ----------------------- User settings -----------------------
dt    = 1e-2;
sim_t = 5.0;

params.m = 1;
params.l = 1;
params.g = 9.81;
params.b = 0.01;
params.I = params.m*params.l^2/3;

params.u_max = 10;
params.u_min = -params.u_max;

params.Kp = 6;
params.Kd = 5;

sigma_linear  = 3;
sigma_concave = 3;
sigma_optimal = 5;

params.sigma_H = sigma_optimal;
params.q_slack = 1e5;   % delta slack penalty
params.q_rho   = 1e2;   % default rho penalty

with_slack_linear  = true;
with_slack_optimal = true;
with_slack_concave = false;

x0 = [pi/4; 0.0];

kmin = 0.1;
kmax = 2.0;
p    = 1.0;
r_vals = 1.0;

eps_ratios = [1e-1, 1e-2, 1e-3];

%% ----------------------- Plot styling -----------------------
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesFontSize',16);
set(groot,'defaultAxesLineWidth',1.2);
set(groot,'defaultAxesFontName','Times New Roman');
set(groot,'defaultTextFontName','Times New Roman');

%% ----------------------- Build CLF -----------------------
params.sigma = sigma_concave;
[P, Vfun, ~] = build_quadratic_clf(params);
c = Vfun(x0);
eps_levels = eps_ratios * c;

%% ----------------------- Simulate linear -----------------------
params.sigma = sigma_linear;

N  = ceil(sim_t/dt) + 1;
nx = numel(x0);

X = zeros(N,nx); t = zeros(N,1);
u = zeros(N-1,1); V = zeros(N-1,1); s = zeros(N-1,1);
rho = nan(N-1,1); solve_t = nan(N-1,1); tu = zeros(N-1,1);

X(1,:) = x0.'; t(1)=0; xk = x0; tk = 0;

for k = 1:N-1
    [uk, sk, Vk, rhok, tk_solve] = clf_qp_linear(xk, params, P, with_slack_linear);

    if isempty(uk), uk = 0; else, uk = uk(1); end
    if isempty(sk), sk = NaN; else, sk = sk(1); end
    if isempty(Vk), Vk = xk.'*P*xk; else, Vk = Vk(1); end
    if isempty(rhok), rhok = nan; else, rhok = rhok(1); end
    if isempty(tk_solve), tk_solve = nan; else, tk_solve = tk_solve(1); end

    uk = min(max(uk, params.u_min), params.u_max);

    tu(k)=tk; u(k)=uk; V(k)=Vk; s(k)=sk; rho(k)=rhok; solve_t(k)=tk_solve;

    ode = @(tt,xx) pendulum_dyn(xx, uk, params);
    [~, xs] = ode45(ode, [0 dt], xk);
    xk = xs(end,:).';
    tk = tk + dt;
    X(k+1,:) = xk.'; t(k+1)=tk;
end

traj_lin.t=t; traj_lin.X=X; traj_lin.tu=tu; traj_lin.u=u; traj_lin.V=V;
traj_lin.slack=s; traj_lin.rho=rho; traj_lin.solve_t=solve_t;

%% ----------------------- Simulate optimal-decay, q_rho = 100 -----------------------
params_opt100 = params;
params_opt100.q_slack = 1e5;
params_opt100.q_rho   = 100;

X = zeros(N,nx); t = zeros(N,1);
u = zeros(N-1,1); V = zeros(N-1,1); s = zeros(N-1,1);
rho = nan(N-1,1); solve_t = nan(N-1,1); tu = zeros(N-1,1);

X(1,:) = x0.'; t(1)=0; xk = x0; tk = 0;

for k = 1:N-1
    [uk, sk, Vk, rhok, tk_solve] = clf_qp_optimal_decay(xk, params_opt100, P, with_slack_optimal);

    if isempty(uk), uk = 0; else, uk = uk(1); end
    if isempty(sk), sk = NaN; else, sk = sk(1); end
    if isempty(Vk), Vk = xk.'*P*xk; else, Vk = Vk(1); end
    if isempty(rhok), rhok = nan; else, rhok = rhok(1); end
    if isempty(tk_solve), tk_solve = nan; else, tk_solve = tk_solve(1); end

    uk = min(max(uk, params_opt100.u_min), params_opt100.u_max);

    tu(k)=tk; u(k)=uk; V(k)=Vk; s(k)=sk; rho(k)=rhok; solve_t(k)=tk_solve;

    ode = @(tt,xx) pendulum_dyn(xx, uk, params_opt100);
    [~, xs] = ode45(ode, [0 dt], xk);
    xk = xs(end,:).';
    tk = tk + dt;
    X(k+1,:) = xk.'; t(k+1)=tk;
end

traj_opt_100.t=t; traj_opt_100.X=X; traj_opt_100.tu=tu; traj_opt_100.u=u; traj_opt_100.V=V;
traj_opt_100.slack=s; traj_opt_100.rho=rho; traj_opt_100.solve_t=solve_t;

%% ----------------------- Simulate optimal-decay, q_rho = 200 -----------------------
params_opt200 = params;
params_opt200.q_slack = 1e5;
params_opt200.q_rho   = 200;

X = zeros(N,nx); t = zeros(N,1);
u = zeros(N-1,1); V = zeros(N-1,1); s = zeros(N-1,1);
rho = nan(N-1,1); solve_t = nan(N-1,1); tu = zeros(N-1,1);

X(1,:) = x0.'; t(1)=0; xk = x0; tk = 0;

for k = 1:N-1
    [uk, sk, Vk, rhok, tk_solve] = clf_qp_optimal_decay(xk, params_opt200, P, with_slack_optimal);

    if isempty(uk), uk = 0; else, uk = uk(1); end
    if isempty(sk), sk = NaN; else, sk = sk(1); end
    if isempty(Vk), Vk = xk.'*P*xk; else, Vk = Vk(1); end
    if isempty(rhok), rhok = nan; else, rhok = rhok(1); end
    if isempty(tk_solve), tk_solve = nan; else, tk_solve = tk_solve(1); end

    uk = min(max(uk, params_opt200.u_min), params_opt200.u_max);

    tu(k)=tk; u(k)=uk; V(k)=Vk; s(k)=sk; rho(k)=rhok; solve_t(k)=tk_solve;

    ode = @(tt,xx) pendulum_dyn(xx, uk, params_opt200);
    [~, xs] = ode45(ode, [0 dt], xk);
    xk = xs(end,:).';
    tk = tk + dt;
    X(k+1,:) = xk.'; t(k+1)=tk;
end

traj_opt_200.t=t; traj_opt_200.X=X; traj_opt_200.tu=tu; traj_opt_200.u=u; traj_opt_200.V=V;
traj_opt_200.slack=s; traj_opt_200.rho=rho; traj_opt_200.solve_t=solve_t;

%% ----------------------- Simulate optimal-decay, q_rho = 300 -----------------------
params_opt300 = params;
params_opt300.q_slack = 1e5;
params_opt300.q_rho   = 300;

X = zeros(N,nx); t = zeros(N,1);
u = zeros(N-1,1); V = zeros(N-1,1); s = zeros(N-1,1);
rho = nan(N-1,1); solve_t = nan(N-1,1); tu = zeros(N-1,1);

X(1,:) = x0.'; t(1)=0; xk = x0; tk = 0;

for k = 1:N-1
    [uk, sk, Vk, rhok, tk_solve] = clf_qp_optimal_decay(xk, params_opt300, P, with_slack_optimal);

    if isempty(uk), uk = 0; else, uk = uk(1); end
    if isempty(sk), sk = NaN; else, sk = sk(1); end
    if isempty(Vk), Vk = xk.'*P*xk; else, Vk = Vk(1); end
    if isempty(rhok), rhok = nan; else, rhok = rhok(1); end
    if isempty(tk_solve), tk_solve = nan; else, tk_solve = tk_solve(1); end

    uk = min(max(uk, params_opt300.u_min), params_opt300.u_max);

    tu(k)=tk; u(k)=uk; V(k)=Vk; s(k)=sk; rho(k)=rhok; solve_t(k)=tk_solve;

    ode = @(tt,xx) pendulum_dyn(xx, uk, params_opt300);
    [~, xs] = ode45(ode, [0 dt], xk);
    xk = xs(end,:).';
    tk = tk + dt;
    X(k+1,:) = xk.'; t(k+1)=tk;
end

traj_opt_300.t=t; traj_opt_300.X=X; traj_opt_300.tu=tu; traj_opt_300.u=u; traj_opt_300.V=V;
traj_opt_300.slack=s; traj_opt_300.rho=rho; traj_opt_300.solve_t=solve_t;

%% ----------------------- Simulate optimal-decay, q_rho = 500 -----------------------
params_opt500 = params;
params_opt500.q_slack = 1e5;
params_opt500.q_rho   = 500;

X = zeros(N,nx); t = zeros(N,1);
u = zeros(N-1,1); V = zeros(N-1,1); s = zeros(N-1,1);
rho = nan(N-1,1); solve_t = nan(N-1,1); tu = zeros(N-1,1);

X(1,:) = x0.'; t(1)=0; xk = x0; tk = 0;

for k = 1:N-1
    [uk, sk, Vk, rhok, tk_solve] = clf_qp_optimal_decay(xk, params_opt500, P, with_slack_optimal);

    if isempty(uk), uk = 0; else, uk = uk(1); end
    if isempty(sk), sk = NaN; else, sk = sk(1); end
    if isempty(Vk), Vk = xk.'*P*xk; else, Vk = Vk(1); end
    if isempty(rhok), rhok = nan; else, rhok = rhok(1); end
    if isempty(tk_solve), tk_solve = nan; else, tk_solve = tk_solve(1); end

    uk = min(max(uk, params_opt500.u_min), params_opt500.u_max);

    tu(k)=tk; u(k)=uk; V(k)=Vk; s(k)=sk; rho(k)=rhok; solve_t(k)=tk_solve;

    ode = @(tt,xx) pendulum_dyn(xx, uk, params_opt500);
    [~, xs] = ode45(ode, [0 dt], xk);
    xk = xs(end,:).';
    tk = tk + dt;
    X(k+1,:) = xk.'; t(k+1)=tk;
end

traj_opt_500.t=t; traj_opt_500.X=X; traj_opt_500.tu=tu; traj_opt_500.u=u; traj_opt_500.V=V;
traj_opt_500.slack=s; traj_opt_500.rho=rho; traj_opt_500.solve_t=solve_t;

%% ----------------------- Simulate concave family -----------------------
params.sigma = sigma_concave;

nr = numel(r_vals);
R = struct('r',[],'traj',[],'label','');

for ir = 1:nr
    r_i = r_vals(ir);

    X = zeros(N,nx); t = zeros(N,1);
    u = zeros(N-1,1); V = zeros(N-1,1); s = zeros(N-1,1);
    rho = nan(N-1,1); solve_t = nan(N-1,1); tu = zeros(N-1,1);

    X(1,:) = x0.'; t(1)=0; xk = x0; tk = 0;

    for k = 1:N-1
        [uk, sk, Vk, rhok, tk_solve] = clf_qp_concave(xk, params, P, c, r_i, kmin, kmax, p, with_slack_concave);

        if isempty(uk), uk = 0; else, uk = uk(1); end
        if isempty(sk), sk = NaN; else, sk = sk(1); end
        if isempty(Vk), Vk = xk.'*P*xk; else, Vk = Vk(1); end
        if isempty(rhok), rhok = nan; else, rhok = rhok(1); end
        if isempty(tk_solve), tk_solve = nan; else, tk_solve = tk_solve(1); end

        uk = min(max(uk, params.u_min), params.u_max);

        tu(k)=tk; u(k)=uk; V(k)=Vk; s(k)=sk; rho(k)=rhok; solve_t(k)=tk_solve;

        ode = @(tt,xx) pendulum_dyn(xx, uk, params);
        [~, xs] = ode45(ode, [0 dt], xk);
        xk = xs(end,:).';
        tk = tk + dt;
        X(k+1,:) = xk.'; t(k+1)=tk;
    end

    traj_i.t=t; traj_i.X=X; traj_i.tu=tu; traj_i.u=u; traj_i.V=V;
    traj_i.slack=s; traj_i.rho=rho; traj_i.solve_t=solve_t;

    R(ir).r = r_i;
    R(ir).traj = traj_i;
    R(ir).label = sprintf('concave, $r=%.1f$, $\\sigma=3$', r_i);
end

%% ----------------------- Quick check -----------------------
disp('max slack (q_rho=100):'); disp(max(traj_opt_100.slack));
disp('max slack (q_rho=200):'); disp(max(traj_opt_200.slack));
disp('max slack (q_rho=300):'); disp(max(traj_opt_300.slack));
disp('max slack (q_rho=500):'); disp(max(traj_opt_500.slack));

disp('max rho   (q_rho=100):'); disp(max(traj_opt_100.rho));
disp('max rho   (q_rho=200):'); disp(max(traj_opt_200.rho));
disp('max rho   (q_rho=300):'); disp(max(traj_opt_300.rho));
disp('max rho   (q_rho=500):'); disp(max(traj_opt_500.rho));

%% ----------------------- Plot style -----------------------
col_lin = [0 0 0];             % linear black
col_con = [1 0 0];             % concave red
col_opt = [0 0 1];   % optimal blue

mk100 = 'o';
mk200 = 's';
mk300 = 'd';
mk500 = '^';

%% ----------------------- Plot (e): V(t) -----------------------
figure('Color','w','Name','Experiment 3: V(t)');
hold on; grid on; box on;
set(gca,'YScale','log');

plot(traj_lin.tu, max(traj_lin.V/c,1e-12), '-', ...
    'Color', col_lin, 'LineWidth', 2.6, ...
    'DisplayName', 'linear, $\sigma=3$');

plot(traj_opt_100.tu, max(traj_opt_100.V/c,1e-12), '-', ...
    'Color', col_opt, 'LineWidth', 2.0, ...
    'Marker', mk100, 'MarkerIndices', 1:20:length(traj_opt_100.tu), ...
    'DisplayName', 'optimal decay, $q_{\rho}=100,\sigma_H=5~~~$');

plot(traj_opt_200.tu, max(traj_opt_200.V/c,1e-12), '--', ...
    'Color', col_opt, 'LineWidth', 2.0, ...
    'Marker', mk200, 'MarkerIndices', 1:20:length(traj_opt_200.tu), ...
    'DisplayName', 'optimal decay, $q_{\rho}=200,\sigma_H=5~~~$');

plot(traj_opt_300.tu, max(traj_opt_300.V/c,1e-12), '-', ...
    'Color', col_opt, 'LineWidth', 2.0, ...
     'MarkerIndices', 1:20:length(traj_opt_300.tu), ...
    'DisplayName', 'optimal decay, $q_{\rho}=300,\sigma_H=5~~~$');

plot(traj_opt_500.tu, max(traj_opt_500.V/c,1e-12), '-.', ...
    'Color', col_opt, 'LineWidth', 2.0, ...
    'Marker', mk500, 'MarkerIndices', 1:20:length(traj_opt_500.tu), ...
    'DisplayName', 'optimal decay, $q_{\rho}=500,\sigma_H=5~~~$');

for ir = 1:nr
    plot(R(ir).traj.tu, max(R(ir).traj.V/c,1e-12), '-', ...
        'Color', col_con, 'LineWidth', 2.2, ...
        'DisplayName', R(ir).label);
end

xlabel('$t$','FontSize',22);
ylabel('$V(x(t))$','FontSize',22);
ylim([1e-4 1]);
yticks([1e-4 1e-3 1e-2 1e-1 1]);
yticklabels({'$10^{-4}V_0$','$10^{-3}V_0$','$10^{-2}V_0$','$10^{-1}V_0$','$V_0$'});
legend('Location','northeast');
legend boxon

%% ----------------------- Plot (f): u versus V -----------------------
figure('Color','w','Name','Experiment 3: u versus V');
hold on; grid off; box on;

% total energies
E_lin = sum(traj_lin.u.^2) * dt;
E_100 = sum(traj_opt_100.u.^2) * dt;
E_200 = sum(traj_opt_200.u.^2) * dt;
E_300 = sum(traj_opt_300.u.^2) * dt;
E_500 = sum(traj_opt_500.u.^2) * dt;

% number of markers on each optimal-decay curve
nmk = 18;

% linear: no marker
plot(traj_lin.V, traj_lin.u, '-', ...
    'Color', col_lin, 'LineWidth', 2.6, ...
    'DisplayName', sprintf('linear, $\\int_0^T u^2\\,dt=%.2f$', E_lin));

% optimal q_rho = 100
Vn = traj_opt_100.V / max(traj_opt_100.V(1),1e-12);
targetV = linspace(1,0,nmk).^2;   % denser at large V, sparser at small V
mk_idx = zeros(1,nmk);
for ii = 1:nmk
    [~, mk_idx(ii)] = min(abs(Vn - targetV(ii)));
end
mk_idx = unique(mk_idx);

plot(traj_opt_100.V, traj_opt_100.u, '-', ...
    'Color', col_opt, 'LineWidth', 2.0, ...
    'Marker', mk100, 'MarkerIndices', mk_idx, ...
    'DisplayName', sprintf('optimal, $q_{\\rho}=100$, $\\int_0^T u^2\\,dt=%.2f$', E_100));

% optimal q_rho = 200
Vn = traj_opt_200.V / max(traj_opt_200.V(1),1e-12);
targetV = linspace(1,0,nmk).^2;
mk_idx = zeros(1,nmk);
for ii = 1:nmk
    [~, mk_idx(ii)] = min(abs(Vn - targetV(ii)));
end
mk_idx = unique(mk_idx);

plot(traj_opt_200.V, traj_opt_200.u, '--', ...
    'Color', col_opt, 'LineWidth', 2.0, ...
    'Marker', mk200, 'MarkerIndices', mk_idx, ...
    'DisplayName', sprintf('optimal, $q_{\\rho}=200$, $\\int_0^T u^2\\,dt=%.2f$', E_200));

% optimal q_rho = 300
Vn = traj_opt_300.V / max(traj_opt_300.V(1),1e-12);
targetV = linspace(1,0,nmk).^2;
mk_idx = zeros(1,nmk);
for ii = 1:nmk
    [~, mk_idx(ii)] = min(abs(Vn - targetV(ii)));
end
mk_idx = unique(mk_idx);

plot(traj_opt_300.V, traj_opt_300.u, '-', ...
    'Color', col_opt, 'LineWidth', 2.0, ...
    'Marker', mk300, 'MarkerIndices', mk_idx, ...
    'DisplayName', sprintf('optimal, $q_{\\rho}=300$, $\\int_0^T u^2\\,dt=%.2f$', E_300));

% optimal q_rho = 500
Vn = traj_opt_500.V / max(traj_opt_500.V(1),1e-12);
targetV = linspace(1,0,nmk).^2;
mk_idx = zeros(1,nmk);
for ii = 1:nmk
    [~, mk_idx(ii)] = min(abs(Vn - targetV(ii)));
end
mk_idx = unique(mk_idx);

plot(traj_opt_500.V, traj_opt_500.u, '-.', ...
    'Color', col_opt, 'LineWidth', 2.0, ...
    'Marker', mk500, 'MarkerIndices', mk_idx, ...
    'DisplayName', sprintf('optimal, $q_{\\rho}=500$, $\\int_0^T u^2\\,dt=%.2f$', E_500));

% concave
for ir = 1:nr
    E_con = sum(R(ir).traj.u.^2) * dt;
    plot(R(ir).traj.V, R(ir).traj.u, '-', ...
        'Color', col_con, 'LineWidth', 2.2, ...
        'DisplayName', sprintf('concave, $r=%.1f$, $\\int_0^T u^2\\,dt=%.2f$', R(ir).r, E_con));
end

xlabel('$V(x)$','FontSize',22);
ylabel('$u$','FontSize',22);
xlim([0 c]);
ylim([0 10]);
legend('Location','best');
legend boxoff

%% ----------------------- Metrics table -----------------------
Controller = {};
EpsRatio = [];
Tconv = [];
Energy = [];
NominalRate = [];
PeakInput = [];
MeanSolve = [];
IntSlack = [];

traj_list = {traj_lin, traj_opt_100, traj_opt_200, traj_opt_300, traj_opt_500};
name_list = {'linear_sigma3', 'optimal_qrho100', 'optimal_qrho200', 'optimal_qrho300', 'optimal_qrho500'};

for ir = 1:numel(R)
    traj_list{end+1} = R(ir).traj;
    name_list{end+1} = sprintf('concave_r%.1f_sigma3', R(ir).r);
end

for ic = 1:numel(traj_list)
    traj = traj_list{ic};
    name = name_list{ic};

    for ie = 1:numel(eps_levels)
        epsV = eps_levels(ie);
        idx = find(traj.V <= epsV, 1, 'first');
        if isempty(idx), idx = numel(traj.V); end

        Tc = traj.tu(idx);
        E  = sum(traj.u(1:idx).^2) * dt;

        V0i = traj.V(1);
        Ve  = traj.V(idx);
        sigma_nom = (1/max(Tc,1e-12)) * log(max(V0i,1e-12)/max(Ve,1e-12));

        st = traj.solve_t(1:idx);
        st = st(~isnan(st));
        if isempty(st), ms = nan; else, ms = mean(st); end

        sl = traj.slack(1:idx);
        sl = sl(~isnan(sl));
        if isempty(sl), sl_int = nan; else, sl_int = sum(sl)*dt; end

        Controller{end+1,1} = name;
        EpsRatio(end+1,1)   = eps_ratios(ie);
        Tconv(end+1,1)      = Tc;
        Energy(end+1,1)     = E;
        NominalRate(end+1,1)= sigma_nom;
        PeakInput(end+1,1)  = max(abs(traj.u(1:idx)));
        MeanSolve(end+1,1)  = ms;
        IntSlack(end+1,1)   = sl_int;
    end
end

Results = table(Controller, EpsRatio, Tconv, Energy, NominalRate, PeakInput, MeanSolve, IntSlack);
disp('--- Experiment 3 metrics ---');
disp(Results);