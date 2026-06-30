function x = Euler_step_nRP2(x, gdot, ddtgdot, De, gamma0, dt,Z)
    % Initialize predictor and corrector variables
    dxdtp = zeros(1, 6);
    dxdtc = zeros(1, 6);
    xp = zeros(1, 6);

    % Define parameters used in calculations
    a1 = -1/9 * x(1) * gdot; % First correction to the stretch
    trQ1 = 2/3 * x(1) * gdot; % Related term with different prefactor

    % Thermodynamic consistency correction
    Lambda = 100;
    beta = -tanh(Lambda * a1);
    beta1 = -Lambda / cosh(Lambda * a1)^2;

    % beta = 1;
    % beta1 = 0;


    % Compute evolution equations for leading order terms (i.e. classic nRP)
    dxdtp(1) = -x(1) + 0 + 6*a1 * (x(1) + beta * (x(1) - 0)) + x(2) * gdot; % xy
    dxdtp(2) = -x(2) + 1 + 6*a1 * (x(2) + beta * (x(2) - 1)) + 0; % yy
    dxdtp(3) = -x(3) + 1 + 6*a1 * (x(3) + beta * (x(3) - 1)) + 2*x(1) * gdot; % xx

    % Compute evolution equations for first correction terms
    a2 = 1/18 * (2/3 * (dxdtp(1) * gdot + x(1) * ddtgdot) - 2*x(4) * gdot + trQ1 - 6*a1 * (trQ1 + beta * trQ1));

    dxdtp(4) = -x(4) + 6*a2 * (x(1) + beta * (x(1) - 0)) + 6*a1 * (x(4) + beta * (a1 * x(1) + x(4) - 0) + beta1 * a2 * (x(1) - 0)) + gdot * x(5); % xy
    dxdtp(5) = -x(5) + 6*a2 * (x(2) + beta * (x(2) - 1)) + 6*a1 * (x(5) + beta * (a1 * x(2) + x(5) - a1) + beta1 * a2 * (x(2) - 1)) + 0; % yy
    dxdtp(6) = -x(6) + 6*a2 * (x(3) + beta * (x(3) - 1)) + 6*a1 * (x(6) + beta * (a1 * x(3) + x(6) - a1) + beta1 * a2 * (x(3) - 1)) + 2 * gdot * x(4); % xx

    % Predictor step
    xp = x + dt * dxdtp;

    % Corrector step
    gdot = gdot + dt * ddtgdot;
    a1 = -1/9 * xp(1) * gdot;
    trQ1 = 2/3 * xp(1) * gdot;

    beta = -tanh(Lambda * a1);
    beta1 = -Lambda / cosh(Lambda * a1)^2;
    % beta = 1;
    % beta1 = 0;

    dxdtc(1) = -xp(1) + 0 + 6*a1 * (xp(1) + beta * (xp(1) - 0)) + xp(2) * gdot; % xy
    dxdtc(2) = -xp(2) + 1 + 6*a1 * (xp(2) + beta * (xp(2) - 1)) + 0; % yy
    dxdtc(3) = -xp(3) + 1 + 6*a1 * (xp(3) + beta * (xp(3) - 1)) + 2*xp(1) * gdot; % xx

    a2 = 1/18 * (2/3 * (dxdtc(1) * gdot + xp(1) * ddtgdot) - 2*xp(4) * gdot + trQ1 - 6*a1 * (trQ1 + beta * trQ1));

    dxdtc(4) = -xp(4) + 6*a2 * (xp(1) + beta * (xp(1) - 0)) + 6*a1 * (xp(4) + beta * (a1 * xp(1) + xp(4) - 0) + beta1 * a2 * (xp(1) - 0)) + gdot * xp(5); % xy
    dxdtc(5) = -xp(5) + 6*a2 * (xp(2) + beta * (xp(2) - 1)) + 6*a1 * (xp(5) + beta * (a1 * xp(2) + xp(5) - a1) + beta1 * a2 * (xp(2) - 1)) + 0; % yy
    dxdtc(6) = -xp(6) + 6*a2 * (xp(3) + beta * (xp(3) - 1)) + 6*a1 * (xp(6) + beta * (a1 * xp(3) + xp(6) - a1) + beta1 * a2 * (xp(3) - 1)) + 2 * gdot * xp(4); % xx

    % Average the predictor and corrector to get the updated x
    x = x + dt * (dxdtp + dxdtc) / 2;



    % Initialize predictor and corrector variables
    % dxdtp = zeros(1, 3);
    % dxdtc = zeros(1, 3);
    % xp = zeros(1, 3);
    % 
    % % Predictor step computations
    % a = sqrt(3 / (2 * x(2) + x(3))); % first correction to the stretch
    % 
    % % Thermodynamic consistency parameters
    % Lambda = 100;
    % fret = 6 * Z * (1 - a);
    % beta = 1; % tanh(Lambda * fret); % Can activate if thermodynamic correction is desired
    % 
    % % Compute predictor evolution equations (classic nRP)
    % dxdtp(1) = -x(1) + 0 - fret * (x(1) + beta * a * (x(1) - 0)) + x(2) * gdot; % xy
    % dxdtp(2) = -x(2) + 1 - fret * (x(2) + beta * a * (x(2) - 1)) + 0;          % yy
    % dxdtp(3) = -x(3) + 1 - fret * (x(3) + beta * a * (x(3) - 1)) + 2 * x(1) * gdot; % xx
    % 
    % % Run the predictor step
    % xp = x + dt * dxdtp;
    % 
    % % Corrector step computations
    % a = sqrt(3 / (2 * xp(2) + xp(3)));
    % 
    % fret = 6 * Z * (1 - a);
    % beta = 1; % tanh(Lambda * fret); % Optional thermodynamic correction
    % gdot = gdot + dt * ddtgdot;
    % 
    % % Compute corrector evolution equations
    % dxdtc(1) = -xp(1) + 0 - fret * (xp(1) + beta * a * (xp(1) - 0)) + xp(2) * gdot; % xy
    % dxdtc(2) = -xp(2) + 1 - fret * (xp(2) + beta * a * (xp(2) - 1)) + 0;            % yy
    % dxdtc(3) = -xp(3) + 1 - fret * (xp(3) + beta * a * (xp(3) - 1)) + 2 * xp(1) * gdot; % xx
    % 
    % % Average predictor and corrector for final update
    % x = x + dt * (dxdtp + dxdtc) / 2;

end
