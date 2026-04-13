
function [u, slack, solve_t] = solve_1d_qp(LfV, LgV, alpha, params, with_slack)
    % Minimize 0.5*u^2 + 0.5*q*slack^2 (if enabled)
    % s.t. LfV + LgV*u <= -alpha + slack
    % and u_min <= u <= u_max

    opts = optimoptions('quadprog','Display','off','ConstraintTolerance',1e-6);

    umin = params.u_min;
    umax = params.u_max;

    t0 = tic;

    if with_slack
        % decision z = [u; d]
        H = diag([1, params.q_slack]);
        f = [0; 0];

        A = [LgV, -1;
             1,    0;
            -1,    0];
        b = [-LfV - alpha;
              umax;
             -umin];

        z = quadprog(H, f, A, b, [], [], [], [], [], opts);
        u = z(1);
        slack = z(2);
    else
        H = 1;
        f = 0;

        A = [LgV;
             1;
            -1];
        b = [-LfV - alpha;
              umax;
             -umin];

        u = quadprog(H, f, A, b, [], [], [], [], [], opts);
        slack = 0;
    end

    solve_t = toc(t0);
end