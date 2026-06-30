% This script runs the LAOS simulation for specific conditions using
% different models (nRP and nRP2). You should run it by section.



%% Run the simulation for Maxwell fluid conditions

% This scripts is the main file to execute the script for performing a stability analysis of a LAOS flow.

clc, clear all, close all
 
% Inputs 
De =2;                    % Deborah number 
gamma0 = 2.5;        % Strain amplitude
Z = 4.5;                   % Entanglement number
% Z =40;                   % Entanglement number
nfm = 8;                  % Number of Fourier modes
ncyc = 20;               % Number of cycles for LAOS simulation
ncyc_sol = 5;           % Number of saved cycles for the analysis
nstep = 1024;          % Number of points per cycle. 


model = @Euler_step_nRP; %(Rolie-Poly model)

% Z = (5:5:20);
for k = 1:length(Z)
    [omega, HTM_eval, HTM_evec] =nRP_LAOS (model, De, gamma0, Z(k), ncyc, ncyc_sol, nfm, nstep);
end



%% Analysis of eigenvector for bulk bands (LTI - Maxwell fluid)

% 
% load test_LTI.mat

[ordered_evals, ordered_evecs]= evals_and_evec_ordering(HTM_eval, HTM_evec); 

clc, clear all, close all
 
% Inputs 
De =1;                    % Deborah number 
gamma0 = 1;        % Strain amplitude
Z = 4;                   % Entanglement number
% Z =40;                   % Entanglement number
nfm = 8;                  % Number of Fourier modes
ncyc = 20;               % Number of cycles for LAOS simulation
ncyc_sol = 5;           % Number of saved cycles for the analysis
nstep = 1024;          % Number of points per cycle. 


model = @Euler_step_nRP; %(Rolie-Poly model)

% Z = (5:5:20);
for k = 1:length(Z)
    [omega, HTM_eval, HTM_evec] =nRP_LAOS (model, De, gamma0, Z(k), ncyc, ncyc_sol, nfm, nstep);
end



%% Analysis of eigenvector for bulk bands (LTI - Maxwell fluid)

% 
% load test_LTI.mat

[ordered_evals, ordered_evecs]= evals_and_evec_ordering(HTM_eval, HTM_evec); 


w= length(omega);
idx = (w-1)/2 +1 ; %center of the frequency band (s =0)

% for i = 1: length(omega)
% 
        V = ordered_evecs {idx}; % Set of nxn eigenvectors associated to the eigenstates in the banda gap
        
        for l = 1: 1: 2*nfm+1
            v  = V (: ,l);    % Localized edge state (eigenvector)
            d = imag(v);            %amplitude of each eigenfunction 
            f = (-nfm:1:nfm);    %Harmonics
            
            figure(10)
            plot(f, d, 'LineStyle','-', 'LineWidth', 1)
            hold on
            xlabel ('Frequency (\omega_n / \omega_o )')
        end 
 % end



idx = (w-1)/2 +1 ; %center of the frequency band (s =0)

% for i = 1: length(omega)
% 
        V = ordered_evecs {idx}; % Set of nxn eigenvectors associated to the eigenstates in the banda gap
        
        for l = 1: 1: 2*nfm+1
            v  = V (: ,l);    % Localized edge state (eigenvector)
            d = abs(v);            %amplitude of each eigenfunction 
            f = (-nfm:1:nfm);    %Harmonics
            
            figure(10)
            plot(f, d, 'LineStyle','-', 'LineWidth', 1)
            hold on
            xlabel ('Frequency (\omega_n / \omega_o )')
        end 
 % end




%% Analysis of eigenvector inside the band gap

clear all, close all

% This matlab file correspond to a simulation with eigenstates inside the
% band gap.

load test.mat

% Eigenvector analysis
evec_analysis (omega, HTM_eval, HTM_evec, nfm, De,nstep)

























%% Eigenvector tracking (before/after EP)
% load test.mat
% 
% [ordered_evals, ordered_evecs ] =evec_track(HTM_eval, HTM_evec);
% 
% 
% 
% 
% [nrows,~] = size(ordered_evals);
% 
% evals_real = real(ordered_evals);
% evals_imag = imag(ordered_evals);
% 
% 
% %Ordered Eigenspectrum
% for n =1:nrows
% plot(evals_real (n,:), evals_imag (n, :), 'o')
% hold on
% % plot(evals_real (n, 1:index), evals_imag (n, 1:index),  'LineStyle', 'none','Marker','x')
% plot(0,0,'+', MarkerSize=10, MarkerFaceColor='r')
% end
%  title('Eigenvalue Spectrum - Ordered  nfm = '+string(nfm) +', De = '+string(De)+ ...
% ', \gamma = '+string(gamma0)+ ', Z = '+string(Z))
% xlabel('Real part')
% ylabel('Imaginary part')
