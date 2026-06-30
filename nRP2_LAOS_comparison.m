function nRP2_LAOS_comparison (fun1, fun2, De, gamma0, Z, ncyc, ncyc_sol, nfm)


%% nRP model (no CCR correction)
[HTM_evals] = nRP_LAOS (fun1, De, gamma0, Z, ncyc, ncyc_sol, nfm);


[nrows, ~] = size(HTM_evals);

ordered_evals= evals_ordering(HTM_evals); 
evals_real = real(ordered_evals);
evals_imag = imag(ordered_evals);


for k =1:nrows
subplot(1,2,1)
plot(evals_real (k, :), evals_imag (k, :), 'o')
hold on
plot(0,0,'+', MarkerSize=10, MarkerFaceColor='r')
% % Adding direction arrows 
% x = evals_real(k,:);
% y = evals_imag(k,:);
% step = 10;  % adjust as needed
% for p = 1:step:length(x)-1
%     dx = x(p+1) - x(p);
%     dy = y(p+1) - y(p);
%         quiver(x(p), y(p), dx, dy, 0, 'k', 'MaxHeadSize', 20);
% end

title('nRP model')
xlabel('Real part')
ylabel('Imaginary part')
end
sgtitle('Eigenvalue Spectrum - Ordered  nfm = '+string(nfm) +', De = '+string(De)+ ...
    ', \gamma = '+string(gamma0)+ ', Z = '+string(Z))



%% nRP2 model (no CCR correction)
[HTM_evals] = nRP2_LAOS (fun2, De, gamma0, Z, ncyc, ncyc_sol, nfm);

[nrows, ~] = size(HTM_evals);

ordered_evals= evals_ordering (HTM_evals); 
evals_real = real(ordered_evals);
evals_imag = imag(ordered_evals);


for k =1:nrows
subplot(1,2,1)
plot(evals_real (k, :), evals_imag (k, :), 'o')
hold on
plot(0,0,'+', MarkerSize=10, MarkerFaceColor='r')
% % Adding direction arrows 
% x = evals_real(k,:);
% y = evals_imag(k,:);
% step = 10;  % adjust as needed
% for p = 1:step:length(x)-1
%     dx = x(p+1) - x(p);
%     dy = y(p+1) - y(p);
%         quiver(x(p), y(p), dx, dy, 0, 'k', 'MaxHeadSize', 20);
% end

title('nRP2 model')
xlabel('Real part')
ylabel('Imaginary part')
end





end 