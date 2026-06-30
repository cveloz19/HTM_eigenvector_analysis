function [] = complex_modulus_new(omega, evals_matrix)

%Note: This works for Nyquist plots with a "single loop"

% This function plots G' and G'' of the material based on the eigenvalue
% calculations. 

%Inputs:
% Frequency range of the baseband components (omega)
% Matrix with the eigenvalues at different frequencies (evals_matrix)
% Cell structure with the set of eigenvectors (evecs_matrix)


%% Ordering Eigenvalues
% ordered_evals= evals_ordering_2 (evals_matrix); 
f = omega(end)*2;

%The eigenvalues are already ordered
ordered_evals =  evals_matrix;
evals_real = real(ordered_evals);
evals_imag = imag(ordered_evals);

%% NOTE: For our case, we need the eigenvalues with negative imaginary part because those are related with positive frequency bands

neg_mask =  any(evals_imag < 0,2);
neg_evals = ordered_evals(neg_mask,:); %Negative imaginary entries

%% Now we need to order the eigenvalues based on the frequency band (harmonic) associated to each case. 
% Each row track the evolution of one eigenvalue. We can sort them using
% the real part of the eigenvalue. The real part decreases at high
% frequency bands.
 plot(real(neg_evals), imag(neg_evals),'o')

norm_matrix = abs(neg_evals);
[nrows, ncolumns] = size (neg_evals);

%Ordering the eigenvalues according to the frequency band
% [~, index] = sortrows(norm_matrix, 1, 'descend');
q = ((ncolumns-1)/2)+1; 
[~, index] = sortrows(norm_matrix, q, 'descend');


sorted_evals = zeros(nrows,ncolumns);


for j = 1: length(index)
    i = index(j);
    sorted_evals (j,:) =neg_evals(i,:);
end 

% Sorted eigenvalues to generated the plots 
% complex viscosity  η = η'-iη''

evals_real = real(sorted_evals);
evals_imag = -1*imag(sorted_evals);

%Plotting complex viscosity
for k = 1:length(index)
    omega_plot = omega+(k-1)*f; 
    figure(20)
    % loglog(omega_plot, evals_real(k,:),'b', 'LineWidth',1, 'LineStyle',':',  'Marker','o', MarkerSize=4)
    plot(omega_plot, evals_real(k,:),'b', 'LineWidth',1, 'LineStyle',':',  'Marker','o', MarkerSize=4)
    hold on
    % loglog(omega_plot, evals_imag(k,:), 'r', 'LineWidth',1, 'LineStyle',':', 'Marker','o', MarkerSize=4)
    plot(omega_plot, evals_imag(k,:), 'r', 'LineWidth',1, 'LineStyle',':', 'Marker','o', MarkerSize=4)
    xlabel('Frequency (\omega)')
    legend('Real(\eta)', '- Img(\eta)^*')
    title('Complex Viscosity \eta' )
end 


%% Compute the complex modulus G =G' + iG''
% G = η*1i(w) where w = s +mf (frequency band)

[n,~] = size(sorted_evals);

for  m = 1:n
    sorted_evals(m,:) = sorted_evals(m,:).*1i.*(omega+(m-1)*f);
end

evals_real = real(sorted_evals);
evals_imag = imag(sorted_evals);

for k = 1:n
    omega_plot = omega+(k-1)*f; 
    figure(21)
    % loglog(omega_plot, evals_real(k,:),'b', 'LineWidth',1, 'LineStyle',':',  'Marker','o', MarkerSize=4)
    plot(omega_plot, evals_real(k,:),'b', 'LineWidth',1, 'LineStyle',':',  'Marker','o', MarkerSize=4)
    hold on
    % loglog(omega_plot, evals_imag(k,:), 'r', 'LineWidth',1, 'LineStyle',':', 'Marker','o', MarkerSize=4)
    plot(omega_plot, evals_imag(k,:), 'r', 'LineWidth',1, 'LineStyle',':', 'Marker','o', MarkerSize=4)
    xlabel('Frequency (\omega)')
    legend('Real (G)', 'Img (G^*)')
    title('Complex Modulus G')
end 


end 
