function limit_cycle = get_LC_solution_nRP2(model, x0,dt, gdot, ddtgdot, De, gamma0, nstep, ncyc, ncyc_sol, Z)
    % This function runs the LAOS flow and extracts specific variables
    % from the last cycle.

    % Initialize
    x = x0;  % Initial conditions
    dt = 2*pi / (2*pi*De * nstep);
    limit_cycle = zeros(6, nstep);
    
    % test=[],
    % test2=[];
    for i = 1:(ncyc * nstep)
        % Perform an Euler step
            % x = Euler_step(x, gdot(mod(i-1, nstep) + 1), ddtgdot(mod(i-1, nstep) + 1), De, gamma0, dt,Z);
            x = model(x, gdot(i), ddtgdot(i), De, gamma0, dt,Z);
            

            % test (i) =gdot(mod(i-1, nstep) + 1);
            % test2(i) =ddtgdot(mod(i-1, nstep) + 1);
             
        % % Store the last cycle
        if i > ((ncyc - ncyc_sol) * nstep)
            idx = i - (ncyc-ncyc_sol) * nstep;
            limit_cycle(:, idx) = x(:)';
            %limit_cycle(2, idx) = x(4);
        end
    end

end
