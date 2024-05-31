function [F,Phi] = objective_smad_trimer(k)

k     = 10.^(k);
nsets = size(k,1);
F     = zeros(nsets,1);
Phi   = zeros(nsets,1); %no penalty
%% Load data to EGFP-Smad2

load("data.mat")

errorbar_lower = (data.errorbar_lower)';
errorbar_upper = (data.errorbar_upper)';
e = errorbar_upper - errorbar_lower;
errorbar_left = e(1:15); errorbar_right = e(16:end);

time = data.time;
tspan1 = time(1:15); tspan2 = time(16:end);

Smad24n = data.Smad24n;
sd = data.Smad24nSD;

% cyt_data = data.cyt_data;
% cyt_data_SB = data.cyt_data_SB;
nuc_data = (data.nuc_data)';
nuc_data_SB = (data.nuc_data_SB)';


% initial conditions
R = 1;
TGFb = 0.066;
Ract = 0;
Rinact = 0;
SB = 1000;
S2c = 60.6;
S2n = 28.5;
pS2c = 0;
pS2n = 0;
S4c = 50.8;
S4n = 50.8;
pS22c = 0;
pS22n = 0;
pS24c = 0;
pS24n = 0;
pS224c = 0;
pS224n = 0; %--nuc. trimer
ES2c = 60.6;
ES2n = 28.5; %nuc. EGFP-SMAD2
EpS2c = 0;
EpS2n = 0; %nuc. EGFP-SMAD2
EpS2cEpS2c = 0;
EpS2nEpS2n = 0; %2*nuc. EGFP-SMAD2
EpS2cpS2c = 0;
EpS2npS2n = 0; %nuc. EGFP-SMAD2
EpS2cS4c = 0;
EpS2nS4n = 0; %nuc. EGFP-SMAD2
EpS2cEpS2cS4c = 0;
EpS2nEpS2nS4n = 0; %2*nuc. EGFP-SMAD2
EpS2cpS2cS4c = 0;
EpS2npS2nS4n = 0; %nuc. EGFP-SMAD2 %--nuc. trimer

y0 = [R TGFb Ract Rinact SB S2c S2n pS2c pS2n S4c S4n pS22c pS22n pS24c pS24n pS224c pS224n ...
    ES2c ES2n EpS2c EpS2n EpS2cEpS2c EpS2nEpS2n EpS2cpS2c EpS2npS2n EpS2cS4c EpS2nS4n EpS2cEpS2cS4c EpS2nEpS2nS4n EpS2cpS2cS4c EpS2npS2nS4n];
for i = 1:nsets

    %     try

    ktgf = k(i,1);
    kphos = k(i,2);
    konSB = k(i,3);
    kon = k(i,4);
    kdephos = k(i,5);
    CIF = k(i,6);
    kin = k(i,7);
    kex = k(i,8);
    koff = k(i,9);
    koffSB = k(i,10);

    ka =[ktgf kphos konSB kon kdephos CIF kin kex koff koffSB];
    % nParams = length(ka);
    %% Calc errors
    %
    % EGFP-Smad2 nuclear concentrations
    %
    sol = ode15s(@Smad_model_trimer,0:0.1:45,y0,[],ka);
    %         if max(sol.x) < 45
    %             error('ode15s failed')
    %         end
    Yn = (deval(sol,tspan1))';

    sol_SB = ode15s(@Smad_model_trimer_SB,45:0.1:150,Yn(end,:),[],ka);
    %         if max(sol_SB.x) < 150
    %             error('ode15s failed')
    %         end
    Yn_SB = (deval(sol_SB,tspan2))';

    Smad2n = (Yn(:,19) + Yn(:,21) +2*Yn(:,23) +Yn(:,25) +Yn(:,27) +2*Yn(:,29)+Yn(:,31))';
    Smad2n_SB = (Yn_SB(:,19) + Yn_SB(:,21) +2*Yn_SB(:,23) +Yn_SB(:,25) +Yn_SB(:,27) +2*Yn_SB(:,29)+Yn_SB(:,31))';

    %
    % EGFP-Smad2 cytoplasmic concentrations
    %
    % [~,Yc] = ode15s(@(t,y)diffun(t,y,ka),tspan1,y0);
    % [~,Yc_SB] = ode15s(@(t,y)diffunSB(t,y,ka),tspan2,Yc(end,:));
    %
    % Smad2c = (Yc(:,7) + Yc(:,9) + Yc(:,12) +Yc(:,14) +2*Yc(:,15))';
    % Smad2c_SB = (Yc_SB(:,7) + Yc_SB(:,9) + Yc_SB(:,12) +Yc_SB(:,14) +2*Yc_SB(:,15))';

    error_all_Smadn = sum(((Smad2n - nuc_data)./errorbar_left).^2,'all')+sum(((Smad2n_SB - nuc_data_SB)./errorbar_right).^2,'all');
    % error_all_Smadc = sum(((Smad2c - cyt_data)./errorbar_left).^2,'all')+sum(((Smad2c_SB - cyt_data_SB)./errorbar_right).^2,'all');
    %
    %     catch ME
    %         disp(ME)
    %         error_all_Smadn = 1e6;
    %     end
    %pSmad24n complex
    %
    % Smad24n_calc = calc_Smad24n(ka);
    % error_Smad_complex = sum(((Smad24n_calc - Smad24n)./sd).^2,'all');
    %% Calc error
    F(i) =  error_all_Smadn;% + error_Smad_complex;%+ sum((ycalc_e(:,:,i) - yvalse).^2) + sum((ycalc_f(:,:,i) - yvalsf).^2);
end

end
