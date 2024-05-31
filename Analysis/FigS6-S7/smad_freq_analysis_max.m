clc
clear
close all

load("Receptor_input.mat")

database_loc = 'G:\My Drive\Research_TAMU\Projects\Smad signaling\Paper\Code\Data\Rmax';
input = total_conc_min;

Files = extractFileLocations(database_loc,'csv');
nFiles = length(Files);


figure
t = tiledlayout('flow');

for i=1:length(Files)

    T = readtable(Files(i));
    T = rmmissing(T);
    T1 = T(T.active_trimer>0.1,:);
    j = str2double(extractBetween(Files(i),"database",".csv"));

    nexttile
    yyaxis left
    fs = 1; %s
    tspan = 0:1:size(input,1)-1; %24 h or 48h
    [p,f] = pspectrum(input(:,j),tspan);
    T = 1./(f*60); %s
    plot(T,pow2db(p),'LineWidth',2)
    yticks(-100:20:20)
    grid on
    
    % ylabel('Power Spectrum (dB)')
    % title('Default Frequency Resolution')

    loc = T1.loc;

    phi = T1.phi;
    NAR = T1.NAR;
    trise = T1.trise;
    trise = trise/60;

    [po,a] = make_pareto(trise,NAR,phi);

    phiP = phi(~a);
    NARP = NAR(~a);
    triseP = trise(~a);
    locP = loc(~a);

    locF = locP(triseP < 6);

    hold on
    yyaxis right
    c = cdfplot(triseP);
    set(c,'LineWidth',2)
    % yticks(0:0.2:1)
    ylabel('')
    title(['Ract = ',num2str(round(mean(input(:,j)),3)),' nM'])

    xticks(0:20:100)
    xtickangle(0)
    xlabel('')
    set(gca,'FontSize',14,'FontName','Arial')
    xlim([0 100])
    title(t,'Rmin')
end

% load("screen_sets_LH.mat");
% load("k_good.mat")
% % k = xb;
% %convert k to min-1 from sec-1; except k(6) which is CIF
% xb([1:5,7:10]) = xb([1:5,7:10])/60;
%
% tspan = 0:1:size(input,1)-1;
% CIF = total_param(:,1);
% PPase = total_param(:,2);
% Smad2 = total_param(:,3);
% Smad4 = total_param(:,4);
%
%
% %%
% close all
% k = 2; i = 1;j = 1;
% % screen = screen_smad_trimer_paramset(input,CIF(locF(k)),PPase(locF(k)),Smad2(locF(k)),Smad4(locF(k)),xb,tspan,i,j);
% screen = screen_smad_trimer_paramset(input,CIF(85),PPase(85),Smad2(85),Smad4(85),xb,tspan,i,j);
%%

function waterplot(s,f,t)
% Waterfall plot of spectrogram
waterfall(f,t,abs(s)'.^2)
set(gca,XDir="reverse",View=[30 50])
xlabel("Frequency (Hz)")
ylabel("Time (s)")
end

%%
function [po,a] = make_pareto(trise,NAR,phi)

% reject_phi = abs(phi) > 1;
% reject_trise = isnan(trise);
% reject_NAR = abs(NAR) < 0;
% keep_data = ~reject_phi & ~reject_NAR & ~reject_trise;
%
% phi_keep = phi(keep_data);
% NAR_keep = NAR(keep_data);
% trise_keep = trise(keep_data)/60;

ep = 0;

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