% This script runs the LAOS simulation for specific conditions using
% different models (nRP and nRP2). You should run it by section.

clc, clear all, close all



% 1. Specify the folder and filename
folderName1 = 'Figures_fig'; %matlab format
folderName2 = 'Figures_jpg'; %jpg format
% 2. Create the folder (if it doesn't exist)
if ~exist(folderName1, 'dir')
    mkdir(folderName1);
end

if ~exist(folderName2, 'dir')
    mkdir(folderName2);
end

% 3. Run the simulation

% Inputs 
De =2;                 % Deborah number 
gamma0 = 1.5;          % Strain amplitude
% Z = 4.5;             % Entanglement number
Z = (0:0.1:2);
nfm = 8;               % Number of Fourier modes
ncyc = 20;             % Number of cycles for LAOS simulation
ncyc_sol = 5;          % Number of saved cycles for the analysis
nstep = 1024;          % Number of points per cycle. 
model = @Euler_step_nRP; % Rolie-Poly model


evals = zeros(2*nfm+1, length(Z));
evecs = cell(length(Z));

frame = 0;
for k = 1:length(Z)

    [omega, HTM_eval, HTM_evec] =nRP_LAOS (model, De, gamma0, Z(k), ncyc, ncyc_sol, nfm, nstep);
 
    evals (:,k) = HTM_eval(:);
    evecs {k} = HTM_evec;


    Nyquist_plot(HTM_eval, nfm, De, Z(k), gamma0, omega);

   % 4. Saving the figures 
    frame = frame +1;
    
    % % Specify file name. 
    % fileName= sprintf('image%06d', frame);
    % 
    % 
    % fullPath1 = fullfile(folderName1, fileName);
    % saveas(gcf, fullPath1, 'fig');
    % 
    % 
    % fullPath2 = fullfile(folderName2, fileName);
    % saveas(gcf, fullPath2, 'jpg');
    close all

end


%% Analysis of eigenvector for energy bands (localization) 
 
[ordered_evals, ordered_evecs]= evals_and_evec_ordering(evals, evecs); 

w = length(omega);
idx = (w-1)/2 +1 ; %center of the frequency band (s = 0)

for i = 1: length(Z)
 
        V = ordered_evecs {idx};    % Set of nxn eigenvectors associated to the eigenstates in the banda gap
        
        for l = 1: 1: 2*nfm+1
            v  = V (: ,l);          % Localized edge state (eigenvector)
            d = abs(v);             % Amplitude of each eigenfunction 
            f = (-nfm:1:nfm);       % Harmonics
            
            figure(1)
            plot(f, d, 'LineStyle','-', 'LineWidth', 1)
            hold on
            xlabel ('Frequency ($\mathrm{\omega_n}$ / $\mathrm{\omega_o}$ )')
        end 
end



%% Analysis of eigenvector inside the band gap
% 
% clear all, close all

% This matlab file correspond to a simulation with eigenstates inside the
% band gap.
% 
% load test.mat
% 
% % Eigenvector analysis
% evec_analysis (omega, HTM_eval, HTM_evec, nfm, De,nstep)

























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
