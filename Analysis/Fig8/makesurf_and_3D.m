clc
clear
close all

f1 = openfig("vary_smad14_at_various_ppase_05-03-2024_Smad4.fig");
view(-43.0806,15.7642)
set(gca,'FontName','Arial','FontSize',24)
% savesurf(f1,'ppase_smadratio')
% f3 = openfig("ppase_smadratio.fig");
exportgraphics(f1,'ppase_smadratio_Smad4.eps','ContentType','vector')