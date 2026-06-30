function [a, b] = jp_fft(f, nfm, Dep)

% T = round(1/Dep);

N = 2^nextpow2(length(f));
% N = length(f); 

% fhat = fft(f,T, 2);

fhat = fft(f,N, 2);
fhat_i = imag(fhat);
fhat_r = real(fhat);

n = floor(nfm); 

% a = fhat_i(2:(n+1))*(-2/N);
a= fhat_i(N:-1:N -(n-1))*2/N;
b = fhat_r(2:(n+1))*2/N;


end 