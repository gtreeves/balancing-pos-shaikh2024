clc
clear
close all

load('screen_sets_LH.mat')
database_loc = '/Users/razeenshaikh/Library/CloudStorage/GoogleDrive-razeen@tamu.edu/My Drive/Research_TAMU/Projects/Smad signaling/trimer_model/19-03-2023/Rmax';

% database_loc = 'G:\My Drive\Research_TAMU\Projects\Smad signaling\trimer_model\19-03-2023\Rmax';
Files = extractFileLocations(database_loc,'csv');
nFiles = length(Files);

i = 1;
T = readtable(Files(i));
T = rmmissing(T);
T1 = T(T.active_trimer>0.1,:);

loc = T1.loc;

LOC = [2142 1698 2884 2758 1020 2036 2062 2850 3020 1966];
% for i=1:length(LOC)
%     text(trise(LOC(i)),NAR(LOC(i)),phi(LOC(i)),num2str(LOC(i)))
% end
% hold on
% scatter3(trise(LOC),NAR(LOC),phi(LOC),25,'red','filled')
% % close all
screen_sets_loc = loc(LOC);

NCPs = total_param(screen_sets_loc,:);
NCPs(5:14,:) = NCPs;
NCPs(1:4,2) = [100 75 50 25];
NCPs(:,4) = [];
% CIF = scatterby(:,1);
% PPase = scatterby(:,2);
% Smad2 = scatterby(:,3);
% Smad4 = scatterby(:,4);
NCP = repelem(NCPs,18,1);

load("Smad24_bins.mat")
Smad24_bins_extended = [0.01 0.1 0.25 0.5 1 1.5 2 3 4 5 7 10 20 30 1e2 5e2 1e3 1e4]';
% Smad24_bins_extended = [Smad24_bins 1e5 0.5e6 1e6 0.5e7 1e7 0.5e8 1e8 0.5e9]';
NCP(:,4) = repmat(Smad24_bins_extended,14,1); %Smad24
NCP(:,5) = NCP(:,3)./NCP(:,4); %Smad4

t = NCP(:,4); NCP(:,4) = NCP(:,5); NCP(:,5) = t; 


% load("NCPs_combinatorial.mat")
NCP(:,1) = ones(252,1)*50; %CIF
NCP(:,3) = ones(252,1)*50; %Smad1
NCP(:,4) = NCP(:,3)./NCP(:,5); %Smad4

save('smad4tot_NCP2_252.mat','NCP')