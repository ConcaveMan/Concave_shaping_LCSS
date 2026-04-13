function xdot = pendulum_dyn(x, u, params)
    psi = x(1);
    w   = x(2);

    m = params.m; l = params.l; g = params.g; b = params.b; I = params.I;

    psi_dot = w;
    w_dot   = (m*g*l/(2*I))*sin(psi) - (b/I)*w - (1/I)*u;

    xdot = [psi_dot; w_dot];
end