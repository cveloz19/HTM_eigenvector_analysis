function x = Euler_step_nRP(x, gdot, ddtgdot, De, gamma0, dt,Z)

    %Initialize predictor and corrector variables
    dxdtp = zeros(1, 3);
    dxdtc = zeros(1, 3);
    xp = zeros(1, 3);

    % Predictor step computations
    a = sqrt(3 / (2 * x(2) + x(3))); % first correction to the stretch

    % Thermodynamic consistency parameters
    Lambda = 100;
    fret = 6 * Z * (1 - a);
    beta = tanh(Lambda * fret); % Can activate if thermodynamic correction is desired

    % Compute predictor evolution equations (classic nRP)
    dxdtp(1) = -x(1) + 0 - fret * (x(1) + beta * a * (x(1) - 0)) + x(2) * gdot; % xy
    dxdtp(2) = -x(2) + 1 - fret * (x(2) + beta * a * (x(2) - 1)) + 0;          % yy
    dxdtp(3) = -x(3) + 1 - fret * (x(3) + beta * a * (x(3) - 1)) + 2 * x(1) * gdot; % xx

    % Run the predictor step
    xp = x + dt * dxdtp;

    % Corrector step computations
    a = sqrt(3 / (2 * xp(2) + xp(3)));

    fret = 6 * Z * (1 - a);
    beta = tanh(Lambda * fret); % Optional thermodynamic correction
    gdot = gdot + dt * ddtgdot;

    % Compute corrector evolution equations
    dxdtc(1) = -xp(1) + 0 - fret * (xp(1) + beta * a * (xp(1) - 0)) + xp(2) * gdot; % xy
    dxdtc(2) = -xp(2) + 1 - fret * (xp(2) + beta * a * (xp(2) - 1)) + 0;            % yy
    dxdtc(3) = -xp(3) + 1 - fret * (xp(3) + beta * a * (xp(3) - 1)) + 2 * xp(1) * gdot; % xx

    % Average predictor and corrector for final update
    x = x + dt * (dxdtp + dxdtc) / 2;

end
