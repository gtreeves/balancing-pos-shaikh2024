%create receptor_input_var_freq4.mat that runs for 100h

clc
clear
close all

load("Receptor_input.mat")

input_max = total_conc_max;
input_min = total_conc_min;

Rmax = mean(input_max,1);
Rmin = mean(input_min,1);

fs = 1; %s
tspan_max = 0:1:size(input_max,1)-1; %24 h or 48h
[p_max,f_max] = pspectrum(input_max,tspan_max);

tspan_min = 0:1:size(input_min,1)-1; %24 h or 48h
[p_min,f_min] = pspectrum(input_min,tspan_min);
% T = 1./(f*60); %s

%%
% figure
% plot(f_min,'*')
% figure
% plot(f_max,'o')

% f_L = mean([f_min(2),f_max(2)]);
f_L = f_max(2);
% f_U = mean([f_min(end-1),f_max(end-1)]);
f_U = 0.45;%f_max(end-5);
f = (logspace(log10(f_L),log10(f_U),100))';
R = [Rmax Rmin]';
tmax = 0:1:100*60*60-1;
F_max = (sin(2*pi().*f*tmax).*std(input_max(:,1)) + Rmax(1))';
save('receptor_input_var_freq4.mat',"F_max")