clc
clear
close all

database_loc = 'G:\My Drive\Research_TAMU\Projects\Smad signaling\Paper\Code\Data\Rmax';
addpath(database_loc)

Files = extractFileLocations(database_loc,'csv');
nFiles = length(Files);
i = 1;
iloc = ['Ract',num2str(i)];

load("screen_sets_LH.mat")
colors = distinguishable_colors(30);
% OrRd = brewermap(10,'OrRd');

T = readtable(Files(i));
% T = readtable("database01.csv");
T = rmmissing(T);

T1 = T(T.active_trimer>0.1,:);
T1v = T1.Variables;
phi = T1.phi;
NAR = T1.NAR;
trise = T1.trise;
loc = T1.loc;
active_trimer = T1.active_trimer;


Smad2_edges = [0 10 50 100 1000];
Smad4_edges = [1 10 50 100 1000];
PPase_edges = [0.01 0.1 1 10 100];
CIF_edges = [0.1 1 10 50 round(max(total_param(loc,1)),0)];

scatterby = total_param(loc,:);

CIF = scatterby(:,1);
PPase = scatterby(:,2);
Smad2 = scatterby(:,3);
Smad4 = scatterby(:,4);

CIF_dis = discretize(CIF,CIF_edges);
PPase_dis = discretize(PPase,PPase_edges);
Smad2_dis = discretize(Smad2,Smad2_edges);
Smad4_dis = discretize(Smad4,Smad4_edges);

for i=1:4
    CIF_group = CIF_dis == i;
    centroid_CIF(i,:) = mean(T1v(CIF_group,2:4),1);

    PPase_group = PPase_dis == i;
    centroid_PPase(i,:) = mean(T1v(PPase_group,2:4),1);

    Smad2_group = Smad2_dis == i;
    centroid_Smad2(i,:) = mean(T1v(Smad2_group,2:4),1);

    Smad4_group = Smad4_dis == i;
    centroid_Smad4(i,:) = mean(T1v(Smad4_group,2:4),1);

end

S = repmat([100,75,50,25],numel(centroid_CIF(:,1)),1);
s = [800;600;400;200];

figure
scatter3(1,0,0,100,"black","filled")

hold on
plot3(centroid_Smad2(:,1),centroid_Smad2(:,2)/60,centroid_Smad2(:,3),LineWidth=2,LineStyle='-.',Color=colors(3,:))
hold on
scatter3(centroid_Smad2(:,1),centroid_Smad2(:,2)/60,centroid_Smad2(:,3),s,colors(3,:),"filled")

hold on
plot3(centroid_Smad4(:,1),centroid_Smad4(:,2)/60,centroid_Smad4(:,3),LineWidth=2,LineStyle='-.',Color=colors(4,:))
hold on
scatter3(centroid_Smad4(:,1),centroid_Smad4(:,2)/60,centroid_Smad4(:,3),s,colors(4,:),"filled")

hold on
plot3(centroid_PPase(:,1),centroid_PPase(:,2)/60,centroid_PPase(:,3),LineWidth=2,LineStyle='-.',Color=colors(2,:))
hold on
scatter3(centroid_PPase(:,1),centroid_PPase(:,2)/60,centroid_PPase(:,3),s,colors(2,:),"filled")

hold on
plot3(centroid_CIF(:,1),centroid_CIF(:,2)/60,centroid_CIF(:,3),LineWidth=2,LineStyle='-.',Color=colors(1,:))
hold on
scatter3(centroid_CIF(:,1),centroid_CIF(:,2)/60,centroid_CIF(:,3),s,colors(1,:),"filled")

hold on
trise_zeros = zeros(4,1);
phi_zer0s = zeros(4,1);
plot3(trise_zeros,centroid_PPase(:,2)/60,centroid_PPase(:,3),LineWidth=2,LineStyle='-.',Color=colors(2,:))
hold on
scatter3(trise_zeros,centroid_PPase(:,2)/60,centroid_PPase(:,3),s,'MarkerFaceColor',colors(2,:))

hold on
plot3(centroid_PPase(:,1),trise_zeros,centroid_PPase(:,3),LineWidth=2,LineStyle='-.',Color=colors(2,:))
hold on
scatter3(centroid_PPase(:,1),trise_zeros,centroid_PPase(:,3),s,'MarkerFaceColor',colors(2,:))

hold on
plot3(centroid_PPase(:,1),centroid_PPase(:,2)/60,trise_zeros,LineWidth=2,LineStyle='-.',Color=colors(2,:))
hold on
scatter3(centroid_PPase(:,1),centroid_PPase(:,2)/60,trise_zeros,s,'MarkerFaceColor',colors(2,:))

legend('optimum point','SMAD2','','SMAD4','','PPase','','CIF','')
xlabel('\phi')
ylabel('trise [min]')
zlabel('NAR')
xticks(0:0.2:1.4)
yticks(0:20:180)
zticks(0:0.2:1.2)
title('Ract = 0.01 nM')
axis square
grid on
set(gca,"FontSize",18,'FontName','Arial')

%%
optimum_point = [1 0 0];
distance_CIF = (centroid_CIF(:,1) - optimum_point(1)).^2 + (centroid_CIF(:,2)./max(centroid_CIF(:,2)) - optimum_point(2)).^2 + (centroid_CIF(:,3) - optimum_point(3)).^2;
distance_PPase = (centroid_PPase(:,1) - optimum_point(1)).^2 + (centroid_PPase(:,2)./max(centroid_PPase(:,2)) - optimum_point(2)).^2 + (centroid_PPase(:,3) - optimum_point(3)).^2;
distance_SMAD2 = (centroid_Smad2(:,1) - optimum_point(1)).^2 + (centroid_Smad2(:,2)./max(centroid_Smad2(:,2)) - optimum_point(2)).^2 + (centroid_Smad2(:,3) - optimum_point(3)).^2;
distance_SMAD4 = (centroid_Smad4(:,1) - optimum_point(1)).^2 + (centroid_Smad4(:,2)./max(centroid_Smad4(:,2)) - optimum_point(2)).^2 + (centroid_Smad4(:,3) - optimum_point(3)).^2;

figure

plot(1:4,distance_SMAD2,'Color',colors(3,:),'LineWidth',2)
hold on
scatter(1:4,distance_SMAD2,s,'MarkerFaceColor',colors(3,:))

hold on
plot(1:4,distance_SMAD4,'Color',colors(4,:),'LineWidth',2)
hold on
scatter(1:4,distance_SMAD4,s,'MarkerFaceColor',colors(4,:))

hold on
plot(1:4,distance_PPase,'Color',colors(2,:),'LineWidth',2)
hold on
scatter(1:4,distance_PPase,s,'MarkerFaceColor',colors(2,:))

hold on
plot(1:4,distance_CIF,'Color',colors(1,:),'LineWidth',2)
hold on
scatter(1:4,distance_CIF,s,'MarkerFaceColor',colors(1,:))

xlim([0.5 5])
xticks(1:4)
xlabel('Clusters')
ylabel('Euclidean Distance of Cluster')
legend('SMAD2','','SMAD4','','PPase','','CIF','')
set(gca,"FontSize",14)
