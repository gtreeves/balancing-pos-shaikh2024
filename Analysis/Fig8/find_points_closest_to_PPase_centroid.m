clc
clear
close all

load('screen_sets_LH.mat')

% database_loc = '/Users/razeenshaikh/Library/CloudStorage/GoogleDrive-razeen@tamu.edu/My Drive/Research_TAMU/Projects/Smad signaling/trimer_model/19-03-2023/Rmax';
database_loc = 'G:\My Drive\Research_TAMU\Projects\Smad signaling\trimer_model\19-03-2023\Rmax';
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

load("data_centroid.mat","centroid_PPase") % phi  trise NAR

for i=1:10
    d1(:,i) = ((trise./max(trise) - centroid_PPase(i,2)./max(centroid_PPase(:,2))).^2 + (NAR - centroid_PPase(i,3)).^2 + (phi - centroid_PPase(i,1)).^2);
    d(:,i) = sqrt(d1(:,i));
    [~,mloc] = min(d(:,i));
    LOC(i) = mloc;
end

hold on
scatter3(trise(LOC),NAR(LOC),phi(LOC))