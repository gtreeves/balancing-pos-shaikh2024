clc
clear
close all

database_loc = 'G:\My Drive\Research_TAMU\Projects\Smad signaling\Paper\Code\Data\Rmax';
addpath(database_loc)

Files = extractFileLocations(database_loc,'csv');

load("screen_sets_LH.mat")
rng('default')
c = distinguishable_colors(50);
c = c([1:8,10,12,14,17:18,21,23:24,26:27,29:30,32:35,37:38,40,44,45:46],:);
colors = c(randperm(size(c,1)),:);

% colors = distinguishable_colors(30);

nFiles = length(Files);
% s = [800;600;400;200];
s = linspace(30,250,10);
f1 = figure;
f2 = figure;
optimum_point = [1 0 0];
for j=1:nFiles



    T = readtable(Files(j));
    % T = readtable("database01.csv");
    T = rmmissing(T);

    T1 = T(T.active_trimer>0.1,:);
    T1v = T1.Variables;
    T1v(:,3) = T1v(:,3)./max(T1v(:,3)); %4/4/2024

    phi = T1.phi;
    NAR = T1.NAR;
    trise = T1.trise;
    loc = T1.loc;
    active_trimer = T1.active_trimer;


%     Smad2_edges = [0 10 50 100 1000];
%     Smad4_edges = [1 10 50 100 1000];
%     PPase_edges = [0.01 0.1 1 10 100];
%     CIF_edges = [0.1 1 10 50 round(max(total_param(loc,1)),0)];

    scatterby = total_param(loc,:);

    CIF = scatterby(:,1);
    PPase = scatterby(:,2);
    Smad2 = scatterby(:,3);
    Smad4 = scatterby(:,4);
    Smad24 = Smad2./Smad4;

    Smad2_edges = logspace(log10(min(Smad2)),log10(max(Smad2)),11);
    Smad4_edges = logspace(log10(min(Smad4)),log10(max(Smad4)),11);
    PPase_edges = logspace(log10(min(PPase)),log10(max(PPase)),11);
    CIF_edges = logspace(log10(min(CIF)),log10(max(CIF)),11);
    Smad24_edges = logspace(log10(min(Smad24)),log10(max(Smad24)),11);

    CIF_dis = discretize(CIF,CIF_edges);
    PPase_dis = discretize(PPase,PPase_edges);
    Smad2_dis = discretize(Smad2,Smad2_edges);
    Smad4_dis = discretize(Smad4,Smad4_edges);
    Smad24_dis = discretize(Smad24,Smad24_edges);
text = 'SMAD24';
    for i=1:10
        Smad2_group = Smad24_dis == i;
        centroid(i,:) = mean(T1v(Smad2_group,2:4),1);

    end
    
    distance = ((centroid(:,1) - optimum_point(1)).^2 + (centroid(:,2)./max(centroid(:,2)) - optimum_point(2)).^2 + (centroid(:,3) - optimum_point(3)).^2).^0.5;
    figure(f1)
    hold on
    plot(1:length(distance),distance,'Color',colors(j,:),'LineWidth',2)
    hold on
    scatter(1:length(distance),distance,s,colors(j,:),"filled")
    hold off
%     close all

    figure(f2)
    hold on
    plot3(centroid(:,2)/3600,centroid(:,3),centroid(:,1),LineWidth=1,LineStyle='-.',Color=colors(j,:))
    hold on
    scatter3(centroid(:,2)/3600,centroid(:,3),centroid(:,1),s,colors(j,:),"filled")
    hold off
%     close all
    clearvars -except text j Files nFiles total_param colors s f1 f2 optimum_point
end
% load("Receptor_input.mat")
% Ract_mean = mean(total_conc_max,1);

figure(f1)
ylim([0.5 1.5])
yticks(0.5:0.25:1.5)
% legend(num2str(bmp'),'Location','eastoutside')
xlabel('Cluster')
ylabel('Euclidean Distance from optimal point')
% xlim([0.5 4.5])
% xticks(1:4)
% colormap('hsv')
title(text)
set(gca,"FontSize",24,'FontName','Arial')
% saveas(gca,'SMAD4_distance.fig')
% close all

figure(f2)

% legend(num2str(round(Ract_mean',4)),'Location','eastoutside')
% hold on 
% scatter3(0,1,0,100,"black","filled")
grid on

xticks(0:5:40)
yticks(0:0.1:1.4)
zticks(0:0.1:3)
xlabel('t_{rise} [h]')
ylabel('NAR')
zlabel('\phi')

% zlabel('NAR')
% colormap('hsv')
title(text)
set(gca,"FontSize",24,'FontName','Arial')
f2.Position = [440 167 808 631];
view(162.7688,16.8331)
% close all
% close(f2)
% saveas(gca,'SMAD4.fig')