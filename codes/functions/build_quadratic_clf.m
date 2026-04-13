function [P, Vfun, gradV] = build_quadratic_clf(params)
    % Uses the same structure as many CLF–QP pendulum examples:
    % A = [0 1; c_bar - Kp/I, -b/I - Kd/I]

    m = params.m; l = params.l; g = params.g; b = params.b; I = params.I;

    c_bar = m*g*l/(2*I);
    b_bar = b/I;

    A = [0 1;
         c_bar - params.Kp/I,  -b_bar - params.Kd/I];

    Q = 3 * eye(2);
    P = lyap(A', Q);

    Vfun  = @(x) x.'*P*x;
    gradV = @(x) 2*P*x;
end
