clc
clear
close all

load('screen_sets_LH.mat')

database_loc = 'G:\My Drive\Research_TAMU\Projects\Smad signaling\Paper\Code\Data\Rmax';

Files = extractFileLocations(database_loc,'csv');
nFiles = length(Files);

i = 1;
T = readtable(Files(i));
T = rmmissing(T);
T1 = T(T.active_trimer>0.1,:);

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
options = fitoptions('cubicinterp');
C = fit([triseP(NARP<0.8),NARP(NARP<0.8)],phiP(NARP<0.8),'cubicinterp');
f1 = figure;
% p1 = plot(C,[triseP(NARP<0.8),NARP(NARP<0.8)],phiP(NARP<0.8)); %plot all
% data points
p1 = plot(C);
p1(1).FaceAlpha = 0.5;
p1(1).LineStyle = 'none';
p1(1).FaceColor = [0.6 0.6 0.6];
% 
xlim([0 180])
ylim([0 2])
zlim([0 2])
xticks(0:20:180)
yticks(0:0.2:2)
zticks(0:0.2:2)
xtickangle(0)
ytickangle(0)
set(gca,'FontSize',24,'FontName','Arial')
xlabel('t_{rise} [min]')
ylabel('NAR')
zlabel('sensitivity coefficient (\phi)')
axis square

load("data_centroid.mat")

c = distinguishable_colors(50);
rng('default')
c = c([5:8,10,17:18,14,21,23:24,26:27,29:30,32:35,37:38,40,44,45:46],:);
colors = c;


plot_this = centroid_PPase; % phi  trise NAR
plot_this_str = 'PPase';
S = repmat([100,75,50,25],numel(centroid_CIF(:,1)),1);
s = linspace(30,250,10);

% figure
hold on
scatter3(0,0,1,100,"black","filled")



LOC = [2142 1698 2884 2758 1020 2036 2062 2850 3020 1966];

screen_sets_loc = loc(LOC);

NCPs = total_param(screen_sets_loc,:);
NCPs(:,4) = [];

NCP = repelem(NCPs,10,1);

load("Smad24_bins.mat")
NCP(:,4) = repmat(Smad24_bins',10,1); %Smad24
NCP(:,5) = NCP(:,3)./NCP(:,4); %Smad4

t = NCP(:,4); NCP(:,4) = NCP(:,5); NCP(:,5) = t; 


TSmadratio = readtable("Fig8A.csv");
TSmadratio = rmmissing(TSmadratio);
TSmadratio = sortrows(TSmadratio,1);

phiSr = TSmadratio.phi;
NARSr = TSmadratio.NAR;
triseSr = TSmadratio.trise;
triseSr = triseSr/60;

% colorsSr = distinguishable_colors(10);
sSr = flip(linspace(10,200,10));
sub = [1 11 21 31 41 51 61 71 81 91];
% figure
% f2 = figure;
for i = [1 3 4 5 6 7 8 10]
figure(f1)
hold on
plot3(triseSr(sub+i-1),NARSr(sub+i-1),phiSr(sub+i-1),LineWidth=2,LineStyle='-.',Color=colors(i,:))
hold on
scatter3(triseSr(sub+i-1),NARSr(sub+i-1),phiSr(sub+i-1),sSr,colors(i,:),"filled","MarkerEdgeColor",'k')

end
view(30.5236,10.2115)

%%
% savesurf(f1,'ppase_smadratio')
% f3 = openfig("ppase_smadratio.fig");
% exportgraphics(f3,'ppase_smadratio.eps','ContentType','vector')
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