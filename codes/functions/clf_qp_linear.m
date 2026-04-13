function [u, slack, V, rho, solve_t] = clf_qp_linear(x, params, P, with_slack)
    % Linear CLF constraint: LfV + LgV u <= -sigma V (+ slack)

    [LfV, LgV, V] = lie_derivatives(x, params, P);
    alpha = params.sigma * V;

    [u, slack, solve_t] = solve_1d_qp(LfV, LgV, alpha, params, with_slack);
    rho = nan;
end