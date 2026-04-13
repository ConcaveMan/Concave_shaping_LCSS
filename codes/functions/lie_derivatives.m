function [LfV, LgV, V] = lie_derivatives(x, params, P)
    V = x.'*P*x;
    gV = 2*P*x;

    f = pendulum_dyn(x, 0, params);
    g = [0; -1/params.I];

    LfV = gV.'*f;
    LgV = gV.'*g;
end
