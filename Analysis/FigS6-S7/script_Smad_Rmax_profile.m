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
Rmax = total_conc_max;
% Rmin = total_conc_min;
%
tspan = 0:1:size(Rmax,1)-1; %24 h or 48h

%unpack parameter sets
param_sets = size(total_param,1);


 %change this of you loop Rnoise
% j = input("Enter which Ract you want to use (numeric: 1 - 10) = ");
% i = input("Enter which parameter set you want to look at (numeric: 1-10000) = ");
% 
% CIF = total_param(i,1);
% PPase = total_param(i,2);
% Smad2 = total_param(i,3);
% Smad4 = total_param(i,4);
j = 1;%input("Enter which Ract you want to use (numeric: 1 - 10) = ");
i = 1;%input("Enter which parameter set you want to look at (numeric: 1-10000) = ");

CIF = 5.7;
PPase = 1;
Smad2 = 178.2;
Smad4 = 101.6;
Rnoise = Rmax(1:length(tspan),j);
disp(['Ract is ', num2str(mean(Rnoise)),' nM '])
disp(['Smad2 is ',num2str(Smad2),' nM ','Smad4 is ',num2str(Smad4),' nM ','CIF is ',num2str(CIF),' nM ',...
    'PPase is',num2str(PPase),' nM '])

screen = screen_smad_trimer_paramset(Rnoise,CIF,PPase,Smad2,Smad4,k,tspan,i,j);
phi = screen.phi_active_trimer;
trise = screen.risetime_active_trimer;
NAR = screen.NAR_active_trimer;
active_trimer = screen.active_trimer;


    
