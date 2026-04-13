function traj = simulate_zoh_linear(x0, dt, T, theta, params, P, with_slack)
    N  = ceil(T/dt) + 1;
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
        [uk, sk, Vk, rhok, tk_solve] = clf_qp_linear(xk, params, P, with_slack);

        uk = min(max(uk, -theta), theta);

        tu(k) = tk;
        u(k)  = uk;
        V(k)  = Vk;
        s(k)  = sk;
        rho(k)= rhok;
        solve_t(k) = tk_solve;

        ode = @(tt,xx) pendulum_dyn(xx, uk, params);
        [~, xs] = ode45(ode, [0 dt], xk);

        xk = xs(end,:).';
        tk = tk + dt;

        X(k+1,:) = xk.';
        t(k+1)   = tk;
    end

    traj.t  = t;
    traj.X  = X;
    traj.tu = tu;
    traj.u  = u;
    traj.V  = V;
    traj.slack = s;
    traj.rho = rho;
    traj.solve_t = solve_t;
end