function [u, slack, V, rho, solve_t] = clf_qp_concave(x, params, P, c, r, kmin, kmax, p, with_slack)
    % Concave CLF constraint: LfV + LgV u <= -sigma s(V)V (+ slack)

    [LfV, LgV, V] = lie_derivatives(x, params, P);
    sV = s_rat(V, c, r, kmin, kmax, p);
    alpha = params.sigma * sV * V;

    [u, slack, solve_t] = solve_1d_qp(LfV, LgV, alpha, params, with_slack);
    rho = nan;
end
