clc
clear
close all

addpath /scratch/user/razeen/trimer/Mats/
addpath /scratch/user/razeen/trimer/Model/
addpath /scratch/user/razeen/trimer/Functions/

%
% Updates:
% Check LHC screen sets,redo
%
load('k_good.mat')
k = xb;
%convert k to min-1 from sec-1; except k(6) which is CIF
k([1:5,7:10]) = k([1:5,7:10])/60;
% k = 0.1*rand(1,10);
%
load('screen_sets_LH.mat')
% load('Receptor_levels_24.mat')
load("Receptor_input.mat")
Rmax = total_conc_max;
% Rmin = total_conc_min;
%
tspan = 0:1:size(Rmax,1)-1; %24 h or 48h

%unpack parameter sets
param_sets = size(total_param,1);


 %change this of you loop Rnoise
CIF = total_param(:,1);
PPase = total_param(:,2);
Smad1 = total_param(:,3);
Smad4 = total_param(:,4);

for j=1%:size(Rmax,2)
    %
    %get avg Rinput from Rmax or Rmin
    %
    Rnoise = Rmax(1:length(tspan),j);
    disp(['Rmax#',num2str(j)])

    %
    %Create a table with Nan's to store screen results
    %
    filename = ['database',num2str(j),'.csv'];
    loc = nan; phi = loc; trise = loc; NAR = loc; active_trimer = loc;
    S2c = nan; S2n = nan; pS2c = nan; pS2n = nan; S4c = nan; S4n = nan; pS22c = nan; pS22n = nan; pS24c = nan; pS24n = nan; pS224c = nan; pS224n  = nan;

    writetable(table(loc,phi,trise,NAR,active_trimer,...
     S2c, S2n, pS2c, pS2n, S4c, S4n, pS22c, pS22n, pS24c, pS24n, pS224c, pS224n),filename);


    parfor i=1:length(total_param)
        loc = i;
%     try
        disp(['parameter#',num2str(i)])
        screen = screen_smad_trimer(Rnoise,CIF(i),PPase(i),Smad1(i),Smad4(i),k,tspan,i,j);
        phi = screen.phi_active_trimer;
        trise = screen.risetime_active_trimer;
        NAR = screen.NAR_active_trimer;
        active_trimer = screen.active_trimer;

        S2c = screen.S2c;
        S2n = screen.S2n;
        pS2c = screen.pS2c;
        pS2n = screen.pS2n;
        S4c = screen.S4c;
        S4n = screen.S4n;
        pS22c = screen.pS22c;
        pS22n = screen.pS22n;
        pS24c = screen.pS24c;
        pS24n = screen.pS24n;
        pS224c = screen.pS224c;
        pS224n = screen.pS224n; %--nuc. trimer

        writetable(table(i,phi,trise,NAR,active_trimer,...
            S2c, S2n, pS2c, pS2n, S4c, S4n, pS22c, pS22n, pS24c, pS24n, pS224c, pS224n), ...
            filename,"WriteVariableNames",false,"WriteRowNames",true,"WriteMode","append");

    end
end
    
    %%
%     clear phi trise tfall
%     for i = 1:length(screen)
%         phi(i,:) = screen(i).phi_pSmad24n;
%         NAR(i,:) = screen(i).NAR_pSmad24n;
%         trise(i,:) = screen(i).risetime_pSmad24n;
%     end
% %     figure
%     % hold on
%     % scatter3(abs(phi),NAR,trise);
%     % 
%     % grid on
%     % xlim([0 1])
%     % xlabel('\phi')
%     % ylabel('NAR')
%     % zlabel('trise (min)')
% %     reject_phi = abs(phi) > 1;
% %     reject_trise = isnan(trise);
% %     reject_NAR = NAR < 0;
% %     keep_data = ~reject_phi & ~reject_NAR & ~reject_trise;
% % 
% %     phi = phi(keep_data);
% %     NAR = NAR(keep_data);
% %     trise = trise(keep_data);
% 
%     f = figure;
%     t = tiledlayout(1,3);
% 
%     nexttile
%     scatter(1-abs(phi),NAR,20,'filled')
%     xlim([0 max(1-abs(phi))])
%     ylim([min(NAR) max(NAR)])
%     xlabel('1-\phi')
%     ylabel('NAR')
%     set(gca,'FontSize',14)
% 
%     nexttile
%     scatter(NAR,trise/60,20,'filled')
%     xlim([min(NAR) max(NAR)])
%     ylim([0 max(trise/60)])
%     xlabel('NAR')
%     ylabel('trise [min]')
%     set(gca,'FontSize',14)
% 
%     nexttile
%     scatter(trise/60,1-abs(phi),20,'filled')
%     xlim([0 max(trise/60)])
%     ylim([0 max(1-abs(phi))])
%     xlabel('trise [min]')
%     ylabel('1-\phi')
%     set(gca,'FontSize',14)
% 
%     title(t,['Ract = ',num2str(round(mean(Rnoise),2)),' nM'])
%     subtitle(t,['#',num2str(j)])
%     f.WindowState = 'fullscreen';
%     saveas(gcf,['Plots/Rmax',num2str(j),'.png'])
% end
% 
