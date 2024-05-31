clc
clear
close all

load("Receptor_input.mat")

tnoise_max = 1:1:86400;
tnoise_min = 1:1:172800;



conversion_factor = 120.44;
%%
f = figure;
f.Position = [331 544 947 303];
% f.WindowState = 'fullscreen';
plot(tnoise_max/3600,round(total_conc_max(:,1)*conversion_factor,0))
title(['Ract = ',num2str(round(mean(total_conc_max(:,1)),2)),' nM'])
yticks(0:2:max(round(total_conc_max(:,1)*conversion_factor,0)))
xticks(0:2:max(tnoise_max/3600))

xlabel('time [h]')
set(gca,'FontSize',24)

exportgraphics(gca,['Rmax_min','.eps'],'ContentType','vector')
%%
f = figure;
f.Position = [331 544 947 303];
% f.WindowState = 'fullscreen';
plot(tnoise_max/3600,round(total_conc_max(:,end)*conversion_factor,0))
title(['Ract = ',num2str(round(mean(total_conc_max(:,end)),2)),' nM'])
yticks(0:50:max(round(total_conc_max(:,end)*conversion_factor,0)))
xticks(0:2:max(tnoise_max/3600))

xlabel('time [h]')
set(gca,'FontSize',24)

exportgraphics(gca,['Rmax_max','.eps'],'ContentType','vector')
%%
f = figure;
f.Position = [331 544 947 303];
% f.WindowState = 'fullscreen';
plot(tnoise_min/3600,round(total_conc_min(:,1)*conversion_factor,0))
title(['Ract = ',num2str(round(mean(total_conc_min(:,1)),3)),' nM'])
yticks(0:1:max(round(total_conc_min(:,1)*conversion_factor,0)))
xticks(0:5:max(tnoise_min/3600))

xlabel('time [h]')
set(gca,'FontSize',24)

exportgraphics(gca,['Rmin_min ','.eps'],'ContentType','vector')
%%
%%
f = figure;
f.Position = [331 544 947 303];
% f.WindowState = 'fullscreen';
plot(tnoise_min/3600,round(total_conc_min(:,end)*conversion_factor,0))
title(['Ract = ',num2str(round(mean(total_conc_min(:,end)),2)),' nM'])
yticks(0:5:max(round(total_conc_min(:,end)*conversion_factor,0)))
xticks(0:5:max(tnoise_min/3600))

xlabel('time [h]')
set(gca,'FontSize',24)

exportgraphics(gca,['Rmin_max ','.eps'],'ContentType','vector')
%%
% ylabel('CDF')
% exportgraphics(gca,['receptor_input','.pdf'],'ContentType','vector','Append',true)
%Ract_labels = ["Ract 1","Ract 2","Ract 3","Ract 4","Ract 5",...
   %  "Ract 6","Ract 7","Ract 8","Ract 9","Ract 10",...
   %  "Ract 11","Ract 12","Ract 13","Ract 14","Ract 15",...
   %  "Ract 16","Ract 17","Ract 18","Ract 19","Ract 20",...
   %  "Ract 21","Ract 22","Ract 23","Ract 24","Ract 25",...
   % "Ract 26","Ract 27","Ract 28","Ract 29","Ract 30"];