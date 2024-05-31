clc
clear
close all

load("Smad24_bins.mat","Smad24_bins");
Smad14ratio = Smad24_bins;

CIF = 50*ones(40,1);
PPase = 0.5*ones(40,1);
Smad1 = 50*ones(40,1);
Smad4 = 50*ones(40,1);

NCPs = [CIF PPase Smad1 Smad4];

%vary CIF
CIFv = logspace(-1,2,10);
NCPs(1:10,1) = CIFv;

%vary PPase
PPasev = logspace(-2,2,10);
NCPs(11:20,2) = PPasev;

%vary Smad1
%Smad14ratio = Smad1/Smad4
Smad1v = Smad14ratio.*50;
NCPs(21:30,3) = Smad1v;

%vary Smad4
%Smad14ratio = Smad1/Smad4
Smad4v = 50./Smad14ratio;
NCPs(31:40,4) = Smad4v;

save('NCPs_freq_NAR.mat','NCPs')

