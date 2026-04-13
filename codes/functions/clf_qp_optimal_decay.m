function [u, slack, V, rho, solve_t] = clf_qp_optimal_decay(x, params, P, with_slack)
    % Optimal-decay CLF-QP:
    % min 0.5*u^2 + 0.5*q_rho*(rho-1)^2 + 0.5*q_slack*delta^2
    % s.t. LfV + LgV*u <= -rho*sigma_H*V + delta

    [LfV, LgV, V] = lie_derivatives(x, params, P);

    opts = optimoptions('quadprog','Display','off','ConstraintTolerance',1e-6);

    umin = params.u_min;
    umax = params.u_max;

    t0 = tic;

    if with_slack
        % z = [u; rho; d]
        H = diag([1, params.q_rho, params.q_slack]);
        f = [0; -params.q_rho; 0];

        A = [LgV, params.sigma_H*V, -1;
             1,   0,                0;
            -1,   0,                0;
             0,   1,                0;
             0,  -1,                0;
             0,   0,               -1];
        b = [-LfV;
              umax;
             -umin;
              1;
              0;
              0];

        z = quadprog(H, f, A, b, [], [], [], [], [], opts);

        u     = z(1);
        rho   = z(2);
        slack = z(3);
    else
        % z = [u; rho]
        H = diag([1, params.q_rho]);
        f = [0; -params.q_rho];

        A = [LgV, params.sigma_H*V;
             1,   0;
            -1,   0;
             0,   1;
             0,  -1];
        b = [-LfV;
              umax;
             -umin;
              1;
              0];

        z = quadprog(H, f, A, b, [], [], [], [], [], opts);

        u     = z(1);
        rho   = z(2);
        slack = 0;
    end

    solve_t = toc(t0);
end