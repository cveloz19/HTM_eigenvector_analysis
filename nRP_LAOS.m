function  [omega, HTM_eval, HTM_evec] = nRP_LAOS(model, De, gamma0, Z, ncyc, ncyc_sol, nfm, nstep)
%This script performs a stability analysis of a LAOS flow using the nRP
%model. It perturbs the system with Fourier modes and examines the response to
%these perturbations.

%% SETUP LAOS FLOW (BASE RESPONSE)

% ncyc                  % Number of cycles for LAOS simulation
% ncyc_sol            % Number of saved cycles for the analysis

f = 2*pi*De;                    % Fundamental frequency
% nstep =1024;                  % N points per cycle 

dt = 2*pi / (f* nstep);     % time step

%Defining the shear-rate or shear strain
gdot = zeros(1, nstep*ncyc);
ddtgdot = zeros(1, nstep*ncyc);

for i = 1:nstep*ncyc
    gdot(i) = 2 * pi * De * gamma0 * sin(f * dt * (i));
    ddtgdot(i) = ((2 * pi * De)^2)* gamma0*cos(f * dt * (i));
    % ddtgdot(i) = (gdot(i) - gdot(i - 1)) / dt;
end
% ddtgdot(end) = (gdot(1) - gdot(end)) / dt;


% Solving the base response of the system
x0 = [0, 1, 1]; % x0 = [0, 1, 1, 0, 0, 0]; Initial condition depends on what ODE system (model) are you solving
base_LC = get_base_LC_solution_nRP(model, x0, dt, gdot, ddtgdot, De, gamma0, nstep, ncyc, ncyc_sol, Z);
% base_sxy = base_LC(1, :) + base_LC(4, :) / Z;
base_sxy = base_LC(1, :);
x0 = base_LC(:, end)';



%Test for compare analytic solution for Maxwell fluid
% tspan =(dt:dt:dt*nstep*ncyc_sol);
% g = gamma0;
% B =g*(f^2)/((1+f^2));
% cxyo = (g*f)/(1+f^2)*(sin(f*tspan)-f*cos(f*tspan))+B*exp(-tspan); % Base flow response
% 
% 
% figure(1)
% plot(tspan, base_sxy)
% hold on

%% PERTURBATION ANALYSIS 


% Dep = De;       %Pertubation Frequency  
pe = 0.03 ;      % Perturbation amplitude

% Defining the contour of our Nyquist plot
% n = 50;
% omega_neg = linspace(-(f-1e-3)/2, -f/(2*n), n+1); %
% omega_pos = linspace(f/(2*n), f/2, n+1);
% omega = [omega_neg 0 omega_pos];

omega = [0];

% omega = -1*[omega_neg 0 omega_pos];


% Initializing the matrix of HTM - evals/evecs
% HTM_eval = zeros (2*nfm+1, length (omega));

HTM_eval = zeros (2*nfm+1, length (omega));
HTM_evec = cell (length(omega), 1);


%% START%%
for iter =1: length (omega)
    
    s =omega(iter);
     
    % Initializing Fourier modes (sin/cos perturbations)
    fourier_modes = zeros(4*nfm+2, nstep*ncyc);                    
    ddt_fourier_modes = zeros(4*nfm+2, nstep*ncyc);             
    
    
    %Setting the fourier modes for the analysis (Perturbations)
    for i = 1: 2*nfm+1  %nfm+1
        for j = 1:nstep*ncyc          
                
                fourier_modes(2*i - 1, j) = pe * sin((s+f*(i-nfm-1)) * dt * (j));
                fourier_modes(2*i, j) = pe * cos((s+f*(i-nfm-1))* dt * (j));

                ddt_fourier_modes(2*i - 1, j) =   pe * (s+f*(i-nfm-1))*cos((s+f*(i-nfm-1)) * dt * (j));
                ddt_fourier_modes(2*i, j) = -pe * (s+f*(i-nfm-1))*sin((s+f*(i-nfm-1))* dt * (j));
            
            
        end
    end
        


% Calculating response of LAOS simulation under perturbations     
[yP] = perturb_response_nRP (model, nfm, x0, dt, gdot, fourier_modes, ddtgdot, ...
    ddt_fourier_modes, De, gamma0, nstep, ncyc,ncyc_sol, Z, base_sxy,s);

 


% Defining period T and time span 
T =2*pi/(f);
to =(nstep*(ncyc_sol-1)+1)*dt;
tspan = (to:dt:to+T-dt);

HTM_i = zeros(2*nfm+1, 2*nfm+1);
% HTM = zeros(2*nfm+1,2*nfm+1);
      
                    index =0;
                        for n = -nfm:nfm     
                            index =index+1;
                            for m = 1:2*nfm+1
                                % n = k+m; %Output harmonic

                                id = -1i*(f*(n)).*tspan;


                                 HTM_i (index, m)  =  sum(yP (m,nstep*(ncyc_sol-1)+1:end).* exp(id)) * dt/pe;
                                 
                            end                     

                       end



                    HTM=HTM_i;

                       % 
                       % Figure of the HTM matrix 
                       if s == 0
                          %matrix = (dSdG_cos-1i*dSdG_sin)/pe;
                          matrix = (HTM);
                          real_num = real(matrix);
                          imag_num = imag(matrix);

                          mag_matrix =sqrt(imag_num.^2+real_num.^2);
                          figure()
                          surf(mag_matrix); %Admittance matrix (Steady-state gain of LTPV system)

                          zlabel('Magnitude')
                          xlabel('Input frequency')
                          ylabel('Output frequency')
                          % s.EdgeColor = 'none';
                          colorbar
                       end 

                % Calculating evals/evec of HTM for each s-value.
                % [evecs, evals] =eig(HTM);
                [evecs, evals] =eigs(HTM, 2*nfm+1);
                lambda = diag(evals);
                HTM_eval (:, iter) = lambda;
                HTM_evec {iter,1} = evecs; 
                % HTM_collection {iter,1} =HTM;

end



%% Nyquist Plot
% num_eval = 2*nfm+1; % Number of evals used to generate the Nyquist plot
Nyquist_plot(HTM_eval, nfm, De, Z, gamma0, omega)

%% Complex modulus
% complex_modulus(omega, HTM_eval)
end




