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

F_max = (sin(2*pi().*f*tspan_max).*std(input_max(:,1)) + Rmax(1))';
%%
figure
plot(tspan_max,F_max)
%%
close all
load('NCPs_freq_NAR.mat','NCPs')
%CIF = 1:10,1
%Ppase = 11;20,2
%Smad1tot = 21:30,3
%Smad4tot = 31:40,4
ncp_legend = NCPs(31:40,4);
ncp_legendS = num2str(ncp_legend,'%.1e');

c = distinguishable_colors(50);
rng('default')
% c = c([5:8,10,17:18,14,21,23:24,26:27,29:30,32:35,37:38,40,44,45:46],:);
colors = c;

T = readtable("database.csv");
T = rmmissing(T);
T1 = sortrows(T,{'locj','loc'},{'ascend','ascend'});

figure
for k=31:1:40
    hold on
    plot(1./(60*f),T1(T1.loc ==k,:).NAR,'LineWidth',2,'Marker','+','Color',colors(k,:))
end
xticks([1e-2 1e-1 1e0 1e1 1e2])
ax = gca;
ax.XScale = 'log';
ax.YScale = 'log';
ylabel('NAR')
xlabel('Period (min)')
title('Smad4_{total}')
legend(ncp_legendS,'Location','best','NumColumns',2)
set(gca,'FontSize',14)
