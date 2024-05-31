clc
clear
close all

TSmadratio = readtable("database1.csv");
TSmadratio = rmmissing(TSmadratio);
TSmadratio = sortrows(TSmadratio,1);

phiSr = TSmadratio.phi;
NARSr = TSmadratio.NAR;
triseSr = TSmadratio.trise;
triseSr = triseSr/60;

colorsSr = distinguishable_colors(10);
sSr = linspace(30,250,10);

figure

for i = 1:10

hold on
plot3(triseSr(i:i+9),NARSr(i:i+9),phiSr(i:i+9),LineWidth=2,LineStyle='-.',Color=colorsSr(i,:))
hold on
scatter3(triseSr(i:i+9),NARSr(i:i+9),phiSr(i:i+9),sSr,colorsSr(i,:),"filled")

end