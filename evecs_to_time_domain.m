function [] = evecs_to_time_domain(ordered_evecs, iter, De, nstep, n, nfm, omega, real_in, img_in)

%% Converting the eigenvector into time domain signal.
%Also it generates a set of figures to track the time domain signal
%associated a each eigenvalue of the inner loop.

% Define time span over one period (or more) 
f = 2*pi*De; 
dt = 1/(De*nstep);   
T = (1/De);
tspan = (0:dt:T*2);
   

% Saving figures. 
% 1. Specify the folder and filename
folderName1 = 'MyFigures1';


% 2. Create the folder (if it doesn't exist)
if ~exist(folderName1, 'dir')
    mkdir(folderName1);
end


frame =0; 

for i = 1: length (iter)    

    sxy = 0;
    for l =  1: n
            V = ordered_evecs{l};
            v = V(:, iter(i));         
            
        for k = 1: 2*nfm+1  
                %Time-varying amplitude (Modulation)
                % s = omega(l);
                % sxy = sxy + v(k).*exp(1i*(f*(k-1-nfm)+s).*tspan);

                %No time-varying amplitude (No modulation)
                sxy = sxy + v(k).*exp(1i*(f*(k-1-nfm)).*tspan);
       end
    
frame = frame +1;       
% Specify file name. 
fileName=sprintf('image%06d.jpg', frame);


fig= figure (1);
fig.Units = 'pixels';
fig.Position = [200 150 900 500];
subplot(1,2,1)
if i  == 1
    plot(real_in(i,l), img_in(i,l), 'bo')
    if l == 1 
        text(real_in(i,l), img_in(i,l)+5e-4,  ' - \omega_o / 2', Color= 'b')
    end
    if l == n 
        text(real_in(i,l), img_in(i,l)-5e-4,  '  \omega_o / 2', Color= 'b')
    end
     
else

    plot(real_in(i,l), img_in(i,l), 'ro')

    if l == 1 
    text(real_in(i,l+15), img_in(i,l)+5e-4,  ' - \omega_o / 2', Color = 'r')
    end
    if l == n 
     text(real_in(i,l-15), img_in(i,l)-5e-4,  '  \omega_o / 2', Color = 'r')
    end


end

hold on
xlabel('Real(\lambda)')
ylabel('Img(\lambda)')
xlim([7e-3 12e-3])
ylim([-4e-3 4e-3])


h = subplot(1,2,2);
cla(h)
plot(tspan/T, real(sxy));
hold on
plot(tspan/T, imag(sxy));
xlabel('t/T (time/period)')
ylabel('\gamma')
legend('Real' , 'Imag')
legend ('Box','off')



% 3. Save the figure as jpg. 
fullPath1 = fullfile(folderName1, fileName);
saveas(gcf, fullPath1, 'jpg');


% Display the directory where the figure is saved. 
disp(['Figure saved to: ', fullPath1]);


    end

end 


end 