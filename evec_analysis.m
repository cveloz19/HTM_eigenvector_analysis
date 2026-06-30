function evec_analysis(omega, HTM_eval, HTM_evec, nfm,De, nstep)

%Inputs

%Outputs

% Eigenvectors are the Fourier modes of the periodic base function.
% Eigenvalues represent the fundamental frequency of the carrier function. 



%% START

n = length(omega);
[ordered_evals, ordered_evecs]= evals_and_evec_ordering(HTM_eval, HTM_evec); 




% Identification of the eigenvalues associated to the inner loop. 
% We need to define the region limited by the outer loop
% This is a temporary solution (Maybe there is a better option)

idx = (n-1)/2 +1; %index of the eigenvalues at s = 0,

evals_real = real(ordered_evals);
evals_imag = imag(ordered_evals);


% These are the eigenvalue coordinates for s = 0 (centroide)
% Note: We can identify what centroids are inside the outer loop
xq = evals_real(:,idx);
yq = evals_imag(:,idx);



% Defining the circular region inside the outer loop.

b = 1e-3;
Rmax = max(xq)-b;
ro =  (Rmax/2)+1.5e-4;


r = linspace(0,2*pi, 2*nfm+1);


xv = ro + (Rmax/2)*cos(r);
yv = 0 + (Rmax/2)*sin(r);


% Inner loops inside the defined region. 
[index, ~] = inpolygon(xq,yq, xv, yv);  % index indicate the eigenvalues associated to the inner loops. 

iter = find(index ==1); % index of the eigenvalue corresponding to Fourier modes (row of the eval matrix)


% Eigenvalues of the inner loops at different values of s.  
real_in = evals_real(index,:);
img_in = evals_imag(index,:);

% Visualization
figure()
plot(xv,yv)
hold on
plot (ro, 0,'k+' )
plot (real_in, img_in,'r+' )
plot (evals_real(~index,:),evals_imag(~index,:), 'b+' )
xlabel('Real(\lambda)')
ylabel('Img(\lambda)')



%% Analysis eigenvector structure
evec_amplitude(ordered_evecs, nfm, iter)




%% Converting the eigenvector into time domain signal.
% evecs_to_time_domain(ordered_evecs, iter, De, nstep, n, nfm, omega, real_in, img_in)



end