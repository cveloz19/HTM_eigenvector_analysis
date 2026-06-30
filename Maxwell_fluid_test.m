function [] = Maxwell_fluid_test(xo, s, f, tspan, base_sxy, sin_sxy, cos_sxy, cxyp_sin, cxyp_cos)


g = 1;
pe = 0.001;      % Perturbation amplitude
m = 1;              %Harmonic


k =(s+f*m); 





%Sine perturbation
C = (g*(f^2)/((1+f^2)));
cxy_ref = (g*f)/(1+f^2)*(sin(f*tspan)-f*cos(f*tspan))+C*exp(-tspan); % Base flow response


cxyo = xo(1);                               %Initial condition (EQ)
B =cxyo + k*pe/((1+k^2))+C;

cxy_p1 = pe/(1+(k)^2)*(sin(k*tspan)-k*cos(k*tspan))+(B-C)*exp(-tspan); % Perturbed flow response
cxy_s =cxy_ref +cxy_p1; %Overall response



% Plotting the results
figure(1)
plot(tspan, cxy_s)
hold on
plot(tspan, sin_sxy)
title('Overall response - Sine')
legend('Analytical', 'Numerical')

figure(2)
plot(tspan, cxy_ref)
hold on
plot(tspan, base_sxy)
title('Base response - Sine')
legend('Analytical', 'Numerical')

figure(3)
plot(tspan, cxy_p1)
hold on
plot(tspan, cxyp_sin)
title('Pertubed response - Sine')
legend('Analytical', 'Numerical')




%Cosine perturbation

phi = cxyo +C - pe/(1+k^2);

cxy_p2 = pe/(1+(k)^2)*(cos(k*tspan)+k*sin(k*tspan))+(phi-C)*exp(-tspan); 
cxy_c = cxy_p2 +cxy_ref;





% Plotting the results
figure(4)
plot(tspan, cxy_c)
hold on
plot(tspan, cos_sxy)
title('Overall response - Cosine')
legend('Analytical', 'Numerical')

figure(5)
plot(tspan, cxy_ref)
hold on
plot(tspan, base_sxy)
title('Base response - Cosine')
legend('Analytical', 'Numerical')

figure(6)
plot(tspan, cxy_p2)
hold on
plot(tspan, cxyp_cos)
title('Pertubed response - Cosine')
legend('Analytical', 'Numerical')









end 