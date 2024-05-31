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
    folderloc = strsplit(Files(i),filesep);
    foldername = char(folderloc(end-1));

    if ~exist(foldername,"dir")
        mkdir(foldername);
    else
        if ~exist([foldername,filesep,'Plots_pretty'],"dir")
            mkdir([foldername,filesep,'Plots_pretty']);
        end
    end

    phi = T2(:,2);
    trise = T2(:,3)/60;
    NAR = T2(:,4);

    [po,a] = make_pareto(trise,NAR,phi); %the result po returns normalized (by max) trise and 1-phi
    % figure
    % figure('pos',[402   544   560   420])
    f1 = figure;
%     f.Position(3:4) = [1200 700];
    plot3(trise,NAR,phi,'Color',[0.6 0.6 0.6],'Marker','.','LineStyle','none')
    hold on
    plot3(trise(~a),NAR(~a),phi(~a),'ko')
    hold on
%     xlim([0 1])
%     ylim([0 1])
%     zlim([0 1])
    grid on
    title(['Ract = ',num2str(round(Ract_mean,4)),' nM'])
    xlabel('t_{rise}')
    ylabel('NAR')
    zlabel('\phi')
    set(gca,'Fontsize',14)
    
%     saveas(gcf,[foldername,filesep,'Plots_pretty',filesep,'Rmax_3D_',num2str(i),'.fig']);
%     saveas(gcf,[foldername,filesep,'Plots_pretty',filesep,'Rmax_3D_',num2str(i),'.svg']);
%     saveas(gcf,[foldername,filesep,'Plots_pretty',filesep,'Rmax_3D_',num2str(i),'.jpeg']);
    % close all
%%
    f = figure;
%     f.Position(3:4) = [1500 500];
%     f.WindowState = "fullscreen";
    t = tiledlayout(1,3);

    nexttile
    scatterhistogram(trise,NAR,'MarkerSize',20,'LineWidth',2,'Color',colors(1,:))
        
    xlim([0 max(trise)]); xlabel('trise [min]')
    ylim([0 max(NAR)]); ylabel('NAR')
%     axis equal
    set(gca,'FontSize',14)

    nexttile
    scatterhistogram(NAR,phi,'MarkerSize',20,'LineWidth',2,'Color',colors(5,:))
    
    xlim([0 max(NAR)]); xlabel('NAR')
    ylim([min(phi) max(phi)]); ylabel('\phi')
    
    set(gca,'FontSize',14)

    nexttile
    scatterhistogram(phi,trise,'MarkerSize',20,'LineWidth',2,'Color',colors(9,:))
    
    xlim([min(phi) max(phi)]); xlabel('\phi')
    ylim([0 max(trise)]); ylabel('trise [min]')
    
    set(gca,'FontSize',14)

    %     title(t,['Ract = ',num2str(round(mean(Rnoise),2)),' nM'])
    title(t,['Ract = ',num2str(round(Ract_mean,4)),' nM'])
    % f.WindowState = 'fullscreen';
%     saveas(gcf,[foldername,filesep,'Plots_pretty',filesep,'Projections',num2str(i),'.fig'])
%     saveas(gcf,[foldername,filesep,'Plots_pretty',filesep,'Projections',num2str(i),'.svg'])
%     saveas(gcf,[foldername,filesep,'Plots_pretty',filesep,'Projections',num2str(i),'.jpeg'])
    % close all
%%
    f = figure;
%     f.Position(3:4) = [1500 500];
%     f.WindowState = "fullscreen";
    t = tiledlayout(1,3);
%     po(:,3) = abs(po(:,3));
    [nota_sort,loc] = sort(~a,'descend');
    triseS = trise(loc);
    NARS = NAR(loc);
    phiS = phi(loc);

    nexttile
    s1 = scatterhistogram(triseS,NARS,'GroupData',(nota_sort));
    s1.Color = {colors(2,:),colors(1,:)};
    s1.MarkerSize = [20 20];
    s1.LineWidth = 2;
    xlim([0 max(trise)]); xlabel('trise [min]')
    ylim([0 max(NAR)]); ylabel('NAR')
    set(gca,'FontSize',14)

    nexttile
    s2 = scatterhistogram(NARS,phiS,'GroupData',(nota_sort));
    s2.Color = {colors(6,:),colors(5,:)};
    s2.MarkerSize = [20 20];
    s2.LineWidth = 2;
    xlim([0 max(NAR)]); xlabel('NAR')
    ylim([min(phi) max(phi)]); ylabel('\phi')
%     xticks(0:0.1:1)
    set(gca,'FontSize',14)

    nexttile
    s3 = scatterhistogram(phiS,triseS,'GroupData',(nota_sort));
    s3.Color = {colors(10,:),colors(9,:)};
    s3.MarkerSize = [20 20];
    s3.LineWidth = 2;
    xlim([min(phi) max(phi)]); xlabel('\phi')
%     xticks(0:0.1:1)
    ylim([0 max(trise)]); ylabel('trise [min]')
    set(gca,'FontSize',14)

    %     title(t,['Ract = ',num2str(round(mean(Rnoise),2)),' nM'])
    title(t,['Ract = ',num2str(round(Ract_mean,4)),' nM'])
    subtitle(t,['points on pareto = ',num2str(sum(~a)),'/',num2str(length(a))])
    set(gca,'FontSize',14)
    % f.WindowState = 'fullscreen';

    saveas(gcf,[foldername,filesep,'Plots_pretty',filesep,'ParetoFront2D_',num2str(i),'.fig'])
    saveas(gcf,[foldername,filesep,'Plots_pretty',filesep,'ParetoFront2D_',num2str(i),'.svg'])
    saveas(gcf,[foldername,filesep,'Plots_pretty',filesep,'ParetoFront2D_',num2str(i),'.jpeg'])
    close all
    clear phi trise NAR po a
end
%% function make pareto

function [po,a] = make_pareto(trise,NAR,phi)

% reject_phi = abs(phi) > 1;
% reject_trise = isnan(trise);
% reject_NAR = abs(NAR) < 0;
% keep_data = ~reject_phi & ~reject_NAR & ~reject_trise;
%
% phi_keep = phi(keep_data);
% NAR_keep = NAR(keep_data);
% trise_keep = trise(keep_data)/60;

ep = 0.1;

po = [trise/max(trise) NAR abs(log10(abs(phi)))];
N = length(po);
a = false(N,1);
tic


for i = 1:N
    % 	dPO = PO - PO(i,:);
    % 	if any(all(dPO < 0,2))
    % 		a(i) = true;
    % 	end

    if ~a(i)
        dpo = po - po(i,:);
        v = all(dpo < -ep,2);
        v2 = all(dpo > ep,2);
        a(v2) = true;
        if any(v)
            a(i) = true;
        end
    end

end

% N_to_plot = 5000;
% Y = rand(N,1);
% b = (trise < 31  & (1-phi) < 0.02);
% a = a | b;
% b = b | Y > N_to_plot/N;

%     save(['Pareto/Mats/Rmax_',num2str(j),'.mat'])



end