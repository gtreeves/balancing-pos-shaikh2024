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
options = fitoptions('cubicinterp','ExtrapolationMethod','none');
C = fit([triseP(NARP<0.8),NARP(NARP<0.8)],phiP(NARP<0.8),'cubicinterp');
f1 = figure;
p1 = plot(C,[triseP(NARP<0.8),NARP(NARP<0.8)],phiP(NARP<0.8));
p1(1).FaceAlpha = 0.5;
p1(1).LineStyle = 'none';
p1(1).FaceColor = [0.6 0.6 0.6];
% 
xlim([0 180])
ylim([0 1])
zlim([0 1])
xticks(0:20:180)
yticks(0:0.2:1.4)
zticks(0:0.2:1.4)
xtickangle(0)
ytickangle(0)
set(gca,'FontSize',20,'FontName','Arial')
xlabel('t_{rise} [min]')
ylabel('NAR')
zlabel('sensitivity coefficient (\phi)')
axis square


select_cube(triseP,NARP,phiP,total_param,locP,3,f1);
p1(2).Visible = "off";
f1.Position = [467 210 1095 743];
view(20.8083,18.4827);
% openfig("surf_fit_2.fig");
% c = drawcuboid;
% wait(c)
% 
% c.Position
%%
% figure
% ch = convhull([triseP,NARP,phiP],Simplify=true);
% TO = triangulation(ch,triseP,NARP,phiP);
% trisurf(TO)
% % plot3(PPase_dis_1(ch(:,1),1),PPase_dis_1(ch(:,2),2),PPase_dis_1(ch(:,3),3))
% % hold on
% % scatter3(PPase_dis_1(:,1),PPase_dis_1(:,2),PPase_dis_1(:,3),'black','filled')
% 
% %%
% figure
% k = boundary(triseP,NARP,phiP);
% trisurf(k,triseP,NARP,phiP)
% %%
% figure
% scatter3(triseP,NARP,phiP)
% hold on
% fill3([20 160 40],[0 0 0.5],[0 0 1],'r','FaceAlpha',0.5)
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

%%

function sets_in_cube = select_cube(X,Y,Z,total_param,loc,n,fHandle)

CIF = log10(total_param(loc,1));
PPase = log10(total_param(loc,2));
Smad1 = log10(total_param(loc,3));
Smad4 = log10(total_param(loc,4));
Smad14 = log10(10.^(Smad1)./10.^(Smad4));

colors = [0 0.4470 0.7410;0.9290 0.6940 0.1250;0.4660 0.6740 0.1880];
for i = 1:n
r1 = drawcuboid('Label',['System ',num2str(i)],'Color',colors(i,:));
wait(r1);
x1 = r1.Position(1);
x2 = x1 + r1.Position(4);
y1 = r1.Position(2);
y2 = y1 + r1.Position(5);
z1 = r1.Position(3);
z2 = z1 + r1.Position(6);

xsets = (X > x1) & (X < x2);
ysets = (Y > y1) & (Y < y2);
zsets = (Z > z1) & (Z < z2);

sets_in_rect = find(xsets.*ysets.*zsets);

f2 = figure;
f2.Position = [0 0 750 500];
t = tiledlayout(2,3);
nexttile
histogram(PPase(sets_in_rect),'BinWidth',0.25,'Normalization','pdf','FaceColor',colors(i,:))
xlim([-2 2])
ylim([0 1])
xticks([-2 -1 0 1 2])
xtickangle(0)
axis square
title('PPase')
set(gca,'FontSize',20,'FontName','Arial')

nexttile
histogram(CIF(sets_in_rect),'BinWidth',0.25,'Normalization','pdf','FaceColor',colors(i,:))
xlim([-1 2])
ylim([0 1])
xticks([-1 0 1 2])
xtickangle(0)
axis square
title('CIF')
set(gca,'FontSize',20,'FontName','Arial')

nexttile(4)
histogram(Smad1(sets_in_rect),'BinWidth',0.25,'Normalization','pdf','FaceColor',colors(i,:))
xlim([0 3])
ylim([0 1])
xticks([-1 0 1 2 3])
xtickangle(0)
axis square
title('Smad1')
set(gca,'FontSize',20,'FontName','Arial')

nexttile(5)
histogram(Smad4(sets_in_rect),'BinWidth',0.25,'Normalization','pdf','FaceColor',colors(i,:))
xlim([-1 3])
ylim([0 1])
xticks([-1 0 1 2 3])
xtickangle(0)
axis square
title('Smad4')
set(gca,'FontSize',20,'FontName','Arial')

nexttile(6)
histogram(Smad14(sets_in_rect),'BinWidth',0.25,'Normalization','pdf','FaceColor',colors(i,:))
xlim([-2 3])
ylim([0 1])
xticks([-3 -2 -1 0 1 2 3])
xtickangle(0)
axis square
title('Smad1/Smad4')
set(gca,'FontSize',20,'FontName','Arial')

title(t,['Non-conserved parameters in system',num2str(i)])
set(gca,'FontSize',20,'FontName','Arial')
figure(fHandle)
end
end