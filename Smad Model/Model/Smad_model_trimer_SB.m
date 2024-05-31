function dcdt = Smad_model_trimer_SB(t,c,k)

R = c(1);
TGFb = c(2);
Ract = c(3);
Rinact = c(4);
SB = c(5);
S2c = c(6);
S2n = c(7);
pS2c = c(8);
pS2n = c(9);
S4c = c(10);
S4n = c(11);
pS22c = c(12);
pS22n = c(13);
pS24c = c(14);
pS24n = c(15);
pS224c = c(16);
pS224n = c(17); %nuc. trimer
ES2c = c(18);
ES2n = c(19); %nuc. EGFP-SMAD2
EpS2c = c(20);
EpS2n = c(21); %nuc. EGFP-SMAD2
EpS2cEpS2c = c(22); 
EpS2nEpS2n = c(23); %2*nuc. EGFP-SMAD2
EpS2cpS2c = c(24); 
EpS2npS2n = c(25); %nuc. EGFP-SMAD2
EpS2cS4c = c(26);
EpS2nS4n = c(27); %nuc. EGFP-SMAD2
EpS2cEpS2cS4c = c(28); 
EpS2nEpS2nS4n = c(29); %2*nuc. EGFP-SMAD2
EpS2cpS2cS4c = c(30);
EpS2npS2nS4n = c(31); %nuc. EGFP-SMAD2 %nuc. trimer


ktgf = k(1);
kphos = k(2);
konSB = k(3);
kon = k(4);
kdephos = k(5);
CIF = k(6);
kin = k(7);
kex = k(8);
koff = k(9);
koffSB = k(10);
PPase = 1;

dR = -ktgf*R*TGFb;
dTGFb = -ktgf*R*TGFb;
dRact = ktgf*R*TGFb - konSB*Ract*SB +koffSB*Rinact;
dRinact = konSB*Ract*SB - koffSB*Rinact;
dSB = -konSB*Ract*SB + koffSB*Rinact;

dS2c = -kphos*S2c*Ract -kin*S2c + kex*S2n;
dS2n = kdephos*pS2n*PPase + kin*S2c -kex*S2n;
dpS2c = kphos*S2c*Ract - 2*kon*pS2c^2 +  2*koff*pS22c - kon*pS2c*S4c +  koff*pS24c - kon*pS2c*pS24c + koff*pS224c - kin*pS2c + kex*pS2n...
    - kon*EpS2c*pS2c + koff*EpS2cpS2c - kon*EpS2cS4c*pS2c + koff*EpS2cpS2cS4c;
dpS2n = -kdephos*pS2n*PPase - 2*kon*pS2n^2+ 2*koff*pS22n - kon*pS2n*S4n + koff*pS24n - kon*pS2n*pS24n + koff*pS224n + kin*pS2c - kex*pS2n...
    - kon*EpS2n*pS2n + koff*EpS2npS2n - kon*EpS2nS4n*pS2n + koff*EpS2npS2nS4n;
dS4c = -kon*pS2c*S4c + koff*pS24c -kon*pS22c*S4c + koff*pS224c -kin*S4c + kin*S4n...
    -kon*EpS2c*S4c + koff*EpS2cS4c -kon*EpS2cEpS2c*S4c + koff*EpS2cEpS2cS4c -kon*EpS2cpS2c*S4c + koff*EpS2cpS2cS4c;
dS4n = -kon*pS2n*S4n + koff*pS24n - kon*pS22n*S4n + koff*pS224n + kin*S4c - kin*S4n ...
    - kon*EpS2n*S4n + koff*EpS2nS4n - kon*EpS2nEpS2n*S4n + koff*EpS2nEpS2nS4n - kon*EpS2npS2n*S4n + koff*EpS2npS2nS4n;
dpS22c = kon*pS2c^2 - koff*pS22c - kon*pS22c*S4c + koff*pS224c - kin*pS22c + kex*pS22n;
dpS22n =kon*pS2n^2- koff*pS22n -kon*pS22n*S4n + koff*pS224n + kin*pS22c -kex*pS22n;
dpS24c = kon*pS2c*S4c - koff*pS24c - kon*pS24c*pS2c + koff*pS224c - kin*pS24c + kex*pS24n;
dpS24n =kon*pS2n*S4n -koff*pS24n -kon*pS24n*pS2n + koff*pS224n + kin*pS24c -kex*pS24n;
dpS224c = kon*pS2c*pS24c - koff*pS224c + kon*pS22c*S4c - koff*pS224c - CIF*kin*pS224c + kex*pS224n;
dpS224n = kon*pS2n*pS24n - koff*pS224n + kon*pS22n*S4n - koff*pS224n + CIF*kin*pS224c -kex*pS224n;

dES2c= - kphos*ES2c*Ract -kin*ES2c + kex*ES2n;
dES2n= kdephos*EpS2n*PPase + kin*ES2c - kex*ES2n;
dEpS2c= kphos*ES2c*Ract - kon*EpS2c*S4c + koff*EpS2cS4c - 2*kon*EpS2c^2 + 2*koff*EpS2cEpS2c - kon*EpS2c*pS2c + koff*EpS2cpS2c - kon*EpS2cS4c*EpS2c + koff*EpS2cEpS2cS4c - kin*EpS2c+ kex*EpS2n;
dEpS2n= - kdephos*EpS2n*PPase - kon*EpS2n*S4n + koff*EpS2nS4n - 2*kon*EpS2n^2 + 2*koff*EpS2nEpS2n - kon*EpS2n*pS2n + koff*EpS2npS2n - kon*EpS2nS4n*EpS2n + koff*EpS2nEpS2nS4n + kin*EpS2c - kex*EpS2n;
dEpS2cEpS2c = kon*EpS2c^2 - koff*EpS2cEpS2c - kon*EpS2cEpS2c*S4c + koff*EpS2cEpS2cS4c - kin*EpS2cEpS2c + kex*EpS2nEpS2n;
dEpS2nEpS2n = kon*EpS2n^2 - koff*EpS2nEpS2n - kon*EpS2nEpS2n*S4n + koff*EpS2nEpS2nS4n + kin*EpS2cEpS2c - kex*EpS2nEpS2n;
dEpS2cpS2c = kon*EpS2c*pS2c -koff*EpS2cpS2c -kon*EpS2cpS2c*S4c + koff*EpS2cpS2cS4c -kin*EpS2cpS2c + kex*EpS2npS2n;
dEpS2npS2n = kon*EpS2n*pS2n -koff*EpS2npS2n -kon*EpS2npS2n*S4n + koff*EpS2npS2nS4n + kin*EpS2cpS2c -kex*EpS2npS2n;
dEpS2cS4c = kon*EpS2c*S4c - koff*EpS2cS4c - kon*EpS2cS4c*EpS2c + koff*EpS2cEpS2cS4c - kon*EpS2cS4c*pS2c + koff*EpS2cpS2cS4c - kin*EpS2cS4c + kex*EpS2nS4n;
dEpS2nS4n = kon*EpS2n*S4n - koff*EpS2nS4n - kon*EpS2nS4n*EpS2n + koff*EpS2nEpS2nS4n - kon*EpS2nS4n*pS2n + koff*EpS2npS2nS4n + kin*EpS2cS4c - kex*EpS2nS4n;
dEpS2cEpS2cS4c = kon*EpS2cEpS2c*S4c - koff*EpS2cEpS2cS4c + kon*EpS2cS4c*EpS2c - koff*EpS2cEpS2cS4c - CIF*kin*EpS2cEpS2cS4c + kex*EpS2nEpS2nS4n;
dEpS2nEpS2nS4n = kon*EpS2nEpS2n*S4n - koff*EpS2nEpS2nS4n + kon*EpS2nS4n*EpS2n - koff*EpS2nEpS2nS4n + CIF*kin*EpS2cEpS2cS4c - kex*EpS2nEpS2nS4n;
dEpS2cpS2cS4c = kon*EpS2cpS2c*S4c - koff*EpS2cpS2cS4c + kon*EpS2cS4c*pS2c - koff*EpS2cpS2cS4c - CIF*kin*EpS2cpS2cS4c + kex*EpS2npS2nS4n;
dEpS2npS2nS4n = kon*EpS2npS2n*S4n - koff*EpS2npS2nS4n + kon*EpS2nS4n*pS2n - koff*EpS2npS2nS4n + CIF*kin*EpS2cpS2cS4c - kex*EpS2npS2nS4n;

dcdt = [dR;dTGFb;dRact;dRinact;dSB;dS2c;dS2n;dpS2c;dpS2n;dS4c;dS4n;dpS22c;dpS22n;dpS24c;dpS24n;dpS224c;dpS224n;...
    dES2c;dES2n;dEpS2c;dEpS2n;dEpS2cEpS2c;dEpS2nEpS2n;dEpS2cpS2c;dEpS2npS2n;dEpS2cS4c;dEpS2nS4n;dEpS2cEpS2cS4c;dEpS2nEpS2nS4n;dEpS2cpS2cS4c;dEpS2npS2nS4n];
end

%%
% dS2c = k(8)*S2n -k(7)*S2c -k(2)*S2c*Ract;
% dGc = k(8)*Gn -k(7)*Gc -k(2)*Gc*Ract;
% dpS2c = k(8)*pS2n -k(7)*pS2c+ k(2)*S2c*Ract -k(4)*pS2c*(S4c+ 2*pS2c+ pGc) + k(9)*(S24c+ 2*S22c+ G2c);
% dpGc = k(8)*pGn - k(7)*pGc + k(2)*Gc*Ract -k(4)*pGc*(S4c+ pS2c+ 2*pGc) + k(9)*(G4c+ G2c+ 2*GGc); 
% dS4c = k(7)*S4n -k(7)*S4c -k(4)*S4c*(pS2c+ pGc) + k(9)*(S24c+ G4c);
% dS24c = k(4)*S4c*pS2c - k(9)*S24c - k(7)*k(6)*S24c;
% dG4c = k(4)*pGc*S4c -k(9)*G4c -k(7)*k(6)*G4c;
% dS22c = k(4)*pS2c*pS2c -k(9)*S22c -k(7)*k(6)*S22c;
% dG2c = k(4)*pGc*pS2c -k(9)*G2c -k(7)*k(6)*G2c;
% dGGc = k(4)*pGc*pGc -k(9)*GGc -k(7)*k(6)*GGc;
% dS2n = k(7)*S2c - k(8)*S2n +  k(5)*pS2n*PPase;
% dGn = k(7)*Gc -k(8)*Gn + k(5)*pGn*PPase;
% dpS2n = k(7)*pS2c - k(8)*pS2n - k(5)*pS2n*PPase - k(4)*pS2n*(S4n+ 2*pS2n+ pGn) + k(9)*(S24n+ 2*S22n+ G2n);
% dpGn = k(7)*pGc - k(8)*pGn - k(5)*pGn*PPase -k(4)*pGn*(S4n+ pS2n+ 2*pGn) + k(9)*(G4n+ G2n+ 2*GGn); 
% dS4n = k(7)*S4c -k(7)*S4n -k(4)*S4n*(pS2n+ pGn) + k(9)*(S24n+ G4n);
% dS24n = k(4)*S4n*pS2n - k(9)*S24n +  k(7)*k(6)*S24c;
% dG4n = k(4)*pGn*S4n -k(9)*G4n + k(7)*k(6)*G4c;
% dS22n = k(4)*pS2n*pS2n -k(9)*S22n + k(7)*k(6)*S22c;
% dG2n = k(4)*pGn*pS2n -k(9)*G2n + k(7)*k(6)*G2c;
% dGGn = k(4)*pGn*pGn -k(9)*GGn + k(7)*k(6)*GGc;