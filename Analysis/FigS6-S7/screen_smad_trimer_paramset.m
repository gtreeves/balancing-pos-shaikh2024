function screen = screen_smad_trimer_paramset(Rnoise,CIF,PPase,Smad1tot,Smad4tot,k,tspan,i,j) 
%intial conditions
%G=EGFP-SMad2
%G4=EGFP-SMad2/Smad4 complex
%G2=EGFP-SMad2/Smad2 complex
%GG=EGFP-SMad2/EGFP-SMad2 complex
% R = 1;
%
%newton's method parameters
conv_tol = 1e-8;
conv_tol_now = 1e-3; %if does not converge; change tolerance to this
nStepsMax = 100;
dxj = 1e-8;
alpha = 1;
%
%

%Smad1 in nucleus and cytoplasm
k(6) = CIF;
Keq = k(7)/k(8);
Smad1totc = Smad1tot/(1+Keq);
Smad1totn = Smad1tot - Smad1totc;
%
%Rdelta is 1% more than Ract, used to calc sensitivity coefficient

Ract = mean(Rnoise);
Rdelta = 1.01*Ract;
%
%initial conditions

S2c = Smad1totc;
S2n = Smad1totn;
pS2c = 0;
pS2n = 0;
S4c = Smad4tot;
S4n = Smad4tot;
pS22c = 0;
pS22n = 0;
pS24c = 0;
pS24n = 0;
pS224c = 0;
pS224n = 0; %--nuc. trimer
% ES2c = Smad1totc;
% ES2n = Smad1totn; %nuc. EGFP-SMAD2
% EpS2c = 0;
% EpS2n = 0; %nuc. EGFP-SMAD2
% EpS2cEpS2c = 0;
% EpS2nEpS2n = 0; %2*nuc. EGFP-SMAD2
% EpS2cpS2c = 0;
% EpS2npS2n = 0; %nuc. EGFP-SMAD2
% EpS2cS4c = 0;
% EpS2nS4n = 0; %nuc. EGFP-SMAD2
% EpS2cEpS2cS4c = 0;
% EpS2nEpS2nS4n = 0; %2*nuc. EGFP-SMAD2
% EpS2cpS2cS4c = 0;
% EpS2npS2nS4n = 0; %nuc. EGFP-SMAD2 %--nuc. trimer
%pack initial conditions
c0 = [S2c S2n pS2c pS2n S4c S4n pS22c pS22n pS24c pS24n pS224c pS224n];
% ...
%     ES2c ES2n EpS2c EpS2n EpS2cEpS2c EpS2nEpS2n EpS2cpS2c EpS2npS2n EpS2cS4c EpS2nS4n EpS2cEpS2cS4c EpS2nEpS2nS4n EpS2cpS2cS4c EpS2npS2nS4n];



%% simulate vanilla delta and calculate steady state
opts = odeset('MaxStep',0.1,'NonNegative',1:12);
[t,c] = ode23tb(@Smad_model_trimer_screen,tspan,c0,opts,k,PPase,Ract);
[t_delta,c_delta] = ode23tb(@Smad_model_trimer_screen,tspan,c0,opts,k,PPase,Rdelta);

active_trimer_dynamic = c(:,12);% + c(:,24) +c(:,26);
active_trimer_dynamic_delta = c_delta(:,12);% + c_delta(:,24) + c_delta(:,26);
c_steady_state = ntnDU(@(c,Ract)Smad_model_trimer_screen(t,c,k,PPase,Ract),c(end,:)',conv_tol,nStepsMax,'',dxj,alpha,Ract);

c_delta_steady_state = ntnDU(@(c,Ract)Smad_model_trimer_screen(t,c,k,PPase,Ract),c_delta(end,:)',conv_tol,nStepsMax,'',dxj,alpha,Rdelta);

% we want to know total Smad24 complex in the nucleus (EGFP + without)
active_trimer = c_steady_state(12);% + c_steady_state(24) + c_steady_state(26);
active_trimer_delta = c_delta_steady_state(12);% + c_delta_steady_state(24) + c_delta_steady_state(26);

%% PO1: sensitivity coefficient
phi = (Ract./active_trimer).*((active_trimer_delta-active_trimer)/(Rdelta-Ract));

%% PO2: NAR
opts = odeset('MaxStep',0.5,'NonNegative',1:12);
[tnoise,cnoise] = ode23tb(@Smad_model_trimer_screen_noise,tspan,c_steady_state,opts,k,PPase,Rnoise);
active_trimer_noise = cnoise(:,12);% + cnoise(:,24) + cnoise(:,26);

%standard noise amplification rate (NAR) defined as the ratio of the coefficient of variation between the output and the input
cv_input = std(Rnoise,0,'all')/mean(Rnoise,'all'); %normalized by N-1; weighting scheme 0
cv_output = std(active_trimer_noise,0,'all')/mean(active_trimer_noise,'all');
NAR = cv_output/cv_input;

%% PO3: rise time without inhibitor
control_info = stepinfo(active_trimer_dynamic,t,active_trimer,'RiseTimeLimits',[0 0.95]);

%% save data
% figure
% 
% tiledlayout(1,2)
% nexttile
% plot(t/60,Smad24n_dynamic,'LineWidth',2)
% xlabel('time [min]')
% ylabel('pSmad24n [nM]')
% set(gca,'FontSize',14)
% 
% nexttile
% plot(t/60,Smad24n_dynamic,'LineWidth',2)
% xlabel('time [min]')
% ylabel('pSmad24n [nM]')
% set(gca,'FontSize',14)

% saveas(gcf,['Profile_',num2str(i),'_.png'])
% if j==1
% 
% f = figure;
% f.WindowState = "fullscreen";
% %     f.Position(3:4) = [1500 1000]; %width height
% tl = tiledlayout(1,3);
% nexttile
figure
plot(t/60,active_trimer_dynamic,'LineWidth',2)
hold on
plot(t_delta/60,active_trimer_dynamic_delta,'LineWidth',2)
tmax = ceil(t(end)/60);
xticks(0:20:tmax)
legend('Ract','Ract_{\delta}')
xlabel('time [h]')
ylabel('pactive [nM]')
title('Dynamic Profile')
set(gca,'FontSize',14)
tmaxnoise = ceil(tnoise(end)/3600);
% nexttile
figure
% yyaxis right
conversion_factor = 120.44;
plot(tnoise/3600,round(Rnoise*conversion_factor,0),'LineWidth',0.5)
xticks(0:2:tmaxnoise)
ylabel('Ract (nM)')
set(gca,'FontSize',14)
% hold on
% yyaxis left
figure
plot(tnoise/3600,active_trimer_noise,'LineWidth',2)
ylabel('active trimer conc. (nM)')
set(gca,'FontSize',14)


xticks(0:2:tmaxnoise)
legend('pSmad24n','Ract')
xlabel('time [h]')
title('w/Noise Input')
set(gca,'FontSize',14)

% nexttile
figure
plot(tnoise(end-3600:end)/3600,active_trimer_noise(end-3600:end),'LineWidth',2)
tminnoise_zoom = round(tnoise(end-3600)/3600,0);
tmaxnoise_zoom = round(tnoise(end)/3600,0);
xticks(tminnoise_zoom:0.2:tmaxnoise_zoom)
xlim([tminnoise_zoom tmaxnoise_zoom])
xlabel('time [h]')
ylabel('pSmad24n [nM]')
title('w/Noise Input zoomed')
set(gca,'FontSize',14)

% title(tl,['Ract = ',num2str(round(Ract,4)),' nM'])
% subtitle(tl,['\phi = ',num2str(round(phi,2)),' ; trise = ',num2str(round(control_info.RiseTime/60,0)),'min; NAR = ',num2str(round(NAR,2))])
set(gca,'FontSize',14)

% if ~exist(['Profiles_',num2str(j)],'dir')
%     mkdir(['Profiles_',num2str(j)])
%     addpath(['Profiles_',num2str(j)])
% end

figure
plot(t/60,active_trimer_dynamic,'LineWidth',2)
hold on
plot(t_delta/60,active_trimer_dynamic_delta,'LineWidth',2)
tmax = ceil(t(end));
xlim([0 20])
% xticks(0:1:tmax)
legend('Ract','Ract_{\delta}')
xlabel('time [min]')
ylabel('pactive [nM]')
title('Dynamic Profile')
set(gca,'FontSize',14)


% saveas(gcf,['Profile_',num2str(i),'_',num2str(j),'.png'])
% clear f
% save(['Profiles_',num2str(j),'/Profile_',num2str(i),'_',num2str(j),'.mat'])
%     close all
% end


screen.phi_active_trimer = phi;
screen.risetime_active_trimer = control_info.RiseTime;    
screen.NAR_active_trimer = NAR;
screen.active_trimer = active_trimer;

screen.S2c = c_steady_state(1);
screen.S2n = c_steady_state(2);
screen.pS2c = c_steady_state(3);
screen.pS2n = c_steady_state(4);
screen.S4c = c_steady_state(5);
screen.S4n = c_steady_state(6);
screen.pS22c = c_steady_state(7);
screen.pS22n = c_steady_state(8);
screen.pS24c = c_steady_state(9);
screen.pS24n = c_steady_state(10);
screen.pS224c = c_steady_state(11);
screen.pS224n = c_steady_state(12); %--nuc. trimer
end