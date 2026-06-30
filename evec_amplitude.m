% This function plots the amplitude of each eigenfunction 
% It could help to identify edge states in the system

%Idea: The eigenstate in the band gap can be localized edge states, where
%the eigenfunction cannot propagate in the bulk states. 


function evec_amplitude(ordered_evecs, nfm, iter)


v = zeros(2*nfm+1,2);


k = length(ordered_evecs);

idx =  (k-1)/2 +1;


    V = ordered_evecs {idx}; % Set of nxn eigenvectors associated to the eigenstates in the banda gap
    
    for l = 1: 2*nfm+1
        
        v  = V (:, l);    % Localized edge state (eigenvector)
        d = abs(v)/norm(v);            %amplitude of each eigenfunction 
        f = (-nfm:1:nfm);    %Harmonics
        
        figure(10)
        plot(f, d', 'LineStyle','-', 'LineWidth', 1.5)
        hold on
        xlabel ('Frequency (\omega_n / \omega_o )')

        
    end








    V = ordered_evecs {idx}; % Set of nxn eigenvectors associated to the eigenstates in the banda gap
    
     
    for l = 1: length(iter)
        
        v  = V (:, iter(l));    % Localized edge state (eigenvector)
        d = abs(v)/norm(v);            %amplitude of each eigenfunction 
        f = (-nfm:1:nfm);    %Harmonics
        
        figure(20)
        plot(f, d', 'LineStyle','-', 'LineWidth', 1.5)
        hold on
        xlabel ('Frequency (\omega_n / \omega_o )')

        
    end











end

