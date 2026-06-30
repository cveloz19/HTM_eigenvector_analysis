function [yP, yo] = perturb_response_nRP2 (model, nfm, x0, dt, gdot, fourier_modes, ddtgdot, ddt_fourier_modes,...
    De, gamma0, nstep, ncyc, ncyc_saved, Z, base_sxy0,base_sxy, s,f)
%+ Description: This script run a perturbation analysis using Fourier modes
     

    yP =   zeros(2*nfm+1, nstep*ncyc_saved);
    yo =   zeros(2*nfm+1, nstep*ncyc_saved);% 



    for i = 1:(2*nfm+1) 

        % sin mode perturbation
            % new_LC = get_LC_solution(x0, dt, gdot + fourier_modes(2*i - 1, :), ...
            %                          ddtgdot + ddt_fourier_modes(2*i - 1, :), De, gamma0, nstep, ncyc, ncyc, Z);

            new_LC = get_LC_solution_nRP2(model, x0, dt, gdot + fourier_modes(2*i - 1, :), ...
                                     ddtgdot + ddt_fourier_modes(2*i - 1, :), De, gamma0, nstep, ncyc, ncyc_saved, Z);

            
            to =((ncyc-ncyc_saved)*nstep+1)*dt;
            tspan = (to:dt:ncyc_saved*nstep*dt+to-dt);
 
            new_sxy0 = new_LC(1, :);
            % new_sxy = new_LC(4, :)/Z;
            new_sxy = new_LC(4, :);

  
            dS_sin0 = new_sxy0 - base_sxy0;
            dS_sin = new_sxy - base_sxy;
           

            % cos mode perturbation

            new_LC = get_LC_solution_nRP2(model, x0, dt, gdot + fourier_modes(2*i, :), ...
                                     ddtgdot + ddt_fourier_modes(2*i, :), De, gamma0, nstep, ncyc, ncyc_saved,Z);



            new_sxy0 = new_LC(1, :);
            % new_sxy = new_LC(4, :)/Z;
            new_sxy = new_LC(4, :);
            
            % cos_sxy = new_sxy;
            dS_cos0 = new_sxy0 - base_sxy0;
            dS_cos = new_sxy - base_sxy;
                    
            
            yMo = dS_cos0 + 1i*dS_sin0;
            yM = dS_cos + 1i*dS_sin;

            yo (i, :) = yMo.*exp(-1i*s.*tspan);
            yP (i, :) = yM.*exp(-1i*s.*tspan);


    end




 % Maxwell_fluid_test(x0, s, f, tspan, base_sxy, sin_sxy, cos_sxy,  dS_sin, dS_cos)



end













