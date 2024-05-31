%
clc
clear
close all
%

load("Receptor_input.mat")
database_loc = 'G:\My Drive\Research_TAMU\Projects\Smad signaling\Paper\Code\Data\Rmax';
colors = brewermap(10,'Paired');
Files = extractFileLocations(database_loc,'csv');
nFiles = length(Files);
% 
scatter_size = 10;
for i = 1:nFiles
% i =1;
    disp(['Running:',num2str(i),'/',num2str(nFiles)])
    Ract_mean = mean(total_conc_max(:,i));

    T = readtable(Files(i));
%     T = readtable("database2.csv");
    T1 = rmmissing(T.Variables);

    [~,uloc] = unique(T1(:,1));
    
    T21 = T1(uloc,:);
    T2 = T21(T21(:,5)>0.1,:);

    phi = T2(:,2);
    trise = T2(:,3)/60;
    NAR = T2(:,4);

    figure;
    scatter(trise,NAR,scatter_size,'MarkerFaceColor',[0.6 0.6 0.6],'MarkerEdgeColor','k','LineWidth',0.1);
        
    xlim([0 180]); xticks(0:20:180); xlabel('trise [min]')
    ylim([0 1.2]); yticks(0:0.2:1.2); ylabel('NAR')
    axis square
    set(gca,'FontSize',18,'FontName','Arial')

    figure;
    scatter(NAR,phi,scatter_size,'MarkerFaceColor',[0.6 0.6 0.6],'MarkerEdgeColor','k','LineWidth',0.1);
        
    xlim([0 1.2]); xticks(0:0.2:1.2); xlabel('NAR')
    ylim([0 1.6]); yticks(0:0.2:1.6); ylabel('\phi')
    axis square
    set(gca,'FontSize',18,'FontName','Arial')

    figure;
    scatter(phi,trise,scatter_size,'MarkerFaceColor',[0.6 0.6 0.6],'MarkerEdgeColor','k','LineWidth',0.1);
        
    xlim([0 1.6]); xticks(0:0.2:1.6); xlabel('\phi')
    ylim([0 180]); yticks(0:20:180); ylabel('trise [min]')
    axis square
    set(gca,'FontSize',18,'FontName','Arial')
end