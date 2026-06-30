function limit_cycle = get_base_LC_solution_nRP2(model, x0,dt, gdot, ddtgdot, De, gamma0, nstep, ncyc, ncyc_sol, Z)
    % This function runs the LAOS flow for several cycles and returns
    % the solution along the last cycle.

    % Initialize
    x = x0;  % Initial conditions
    dt = 2*pi / (2*pi*De * nstep);
    % limit_cycle = zeros(6, nstep);
    limit_cycle = zeros(6, nstep);
    
    for i = 1:(ncyc * nstep)
        % Perform an Euler step
        % x = Euler_step(x, gdot(mod(i-1, nstep) + 1), ddtgdot(mod(i-1, nstep) + 1), De, gamma0, dt,Z);
        x = model(x, gdot(i), ddtgdot(i), De, gamma0, dt,Z);
        
        % Store the last cycle (ncyc_sol = 1) 
         if i > ((ncyc-ncyc_sol) * nstep)
            idx = i - (ncyc-ncyc_sol) * nstep;
            limit_cycle(:, idx) = x(:)';
        end
    end

% 
% figure(1)
% plot(gdot(1: ncyc*nstep))
% 
% figure(2)
% plot(ddtgdot(1: ncyc*nstep))

end

