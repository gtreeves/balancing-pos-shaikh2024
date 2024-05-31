clc
clear 
close all


% this code for genarating the .txt file for run parameter sweep on brown
% system

%% input-----------------------

paraset=10000;
parameter_n = 4; % set up for total parameter number
sample_density = 10000; % set up for how many samples 
sample_n = 1; % szl_prameter_sampling number 

%not hypercube
screen_vars = 4;
%% -------------------------------
% latin hypercube sampling
rng('shuffle');
X = lhsdesign(paraset,parameter_n);

CIF = logspace(-1,2,sample_density);
PPase = logspace(-2,2,sample_density);
Smad2 = logspace(-1,3,sample_density);
Smad4 = logspace(-1,3,sample_density);


total_param_sets = [];
for i=1:paraset
% set the putput path name
sample = round(X(i,:).*sample_density);
% check the sample number is 0 if 0 set to first sample
foo = (sample == 0);
sample(foo)=1;
% parameters.Ractset = Ract(sample(1));
parameters.CIFset = CIF(sample(1));
parameters.PPaseset = PPase(sample(2));
parameters.Smad2set = Smad2(sample(3)); 
parameters.Smad4set = Smad4(sample(4)); 
% parameters.kon = kon(sample(4)); 

szl_prameter_sampling = 1;

param = [parameters.CIFset,parameters.PPaseset,parameters.Smad2set,parameters.Smad4set];
% param = [CIF;PPase;Smad1;Smad4];

total_param_sets = [total_param_sets;param];




end
%making total_param from total_param_sets
%initialize
% total_param = zeros(paraset,screen_vars);
% %CIF remains the same
% total_param(:,1) = total_param_sets(:,1);
% %PPase remains the same
% total_param(:,2) = total_param_sets(:,2);
% %Smad2 is created using Smad4 and Smad2/Smad4 ratio in (10e-1,10e1)
% % total_param = [total_param total_param(:,3)];
% factorlog = -1 + (1 + 1)*rand(1,paraset);
% factor = 10.^factorlog;
% total_param(:,3)=factor'.*total_param_sets(:,3);
% %Smad 4 remains the same
% total_param(:,4) = total_param_sets(:,3);


total_param = total_param_sets;
figure
parallelplot(log10(total_param))
save('screen_sets_LH.mat','total_param')

%% 
% cdfplot(R)
% hold on
% cdfplot(PPase)
% hold on
% cdfplot(Smad2)
% hold on
% cdfplot(Smad4)