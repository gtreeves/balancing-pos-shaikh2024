clc
clear
close all

load('screen_sets_LH.mat')
PPase = total_param(:,2);

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


f1 = figure;
% scatter(PPase(locP),triseP)
%{
% options = fitoptions('cubicinterp');
C = fit(PPase(locP),triseP,);
hold on
% p1 = plot(C,[triseP(NARP<0.8),NARP(NARP<0.8)],phiP(NARP<0.8)); %plot all
% data points
p1 = plot(C);
%}
DT = delaunayTriangulation(PPase(locP),triseP);
% Compute the convex hull.

C = convexHull(DT);
% Plot the triangulation and highlight the convex hull in red.
hold on
% plot(DT.Points(:,1),DT.Points(:,2),'.','MarkerSize',10)
hold on
plot(DT.Points(C,1),DT.Points(C,2),'k',LineWidth=2) 
xlim([1e-3 1e3])
xticks([1e-3 1e-2 1e-1 1e0 1e1 1e2 1e3])
ax = gca;
ax.XScale = 'log';
% ax.YScale = 'log';
ylabel('t_{rise} [min]')
xlabel('PPase [nM]')
% title('Smad4_{total}')
% legend(ncp_legendS,'Location','best','NumColumns',2)
set(gca,'FontSize',14,'FontName','Arial')
exportgraphics(f1,'DLconvexHull.eps','ContentType','vector')
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