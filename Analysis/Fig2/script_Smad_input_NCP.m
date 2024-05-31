clc
clear
close all

%
load('k_good.mat')
k = xb;
%convert k to min-1 from sec-1; except k(6) which is CIF
k([1:5,7:10]) = k([1:5,7:10])/60;
% k = 0.1*rand(1,10);
%
load('screen_sets_LH.mat')
% load('Receptor_levels_24.mat')
load("Receptor_input.mat")

useMax = input('Enter 1 to use Rmax inputs, 0 otherwise');

if useMax
    Ract = total_conc_max;
else
    Ract = total_conc_min;
end

%
tspan = 0:1:size(Ract,1)-1; %24 h or 48h

%unpack parameter sets
param_sets = size(total_param,1);
j = input('Choose an Ract input (numeric:1-10) = ');
CIF = input('Input a value for CIF (nM) = ');
PPase = input('Input a value for Ppase (nM) = ');
Smad2 = input('Input a value for Smad2 (nM) = ');
Smad4 = input('Input a value for Smad4 (nM) = ');

Rnoise = Ract(1:length(tspan),j);
disp(['Ract is ', num2str(mean(Rnoise)),' nM '])
disp(['Smad2 is ',num2str(Smad2),' nM ','Smad4 is ',num2str(Smad4),' nM ','CIF is ',num2str(CIF),' nM ',...
    'PPase is ',num2str(PPase),' nM '])
i = abs(randi(5000));
screen = screen_smad_trimer_paramset(Rnoise,CIF,PPase,Smad2,Smad4,k,tspan,i,j);
phi = screen.phi_active_trimer;
trise = screen.risetime_active_trimer;
NAR = screen.NAR_active_trimer;
active_trimer = screen.active_trimer;


    
