function s = s_rat(V, c, r, kmin, kmax, p)
    % Rational shape with s(c)=r
    c = max(c, eps);
    ell = (r - kmin) / max(kmax - r, 1e-12);
    z = (V./c).^p;
    s = (kmin*z + kmax*ell) ./ (z + ell);
end