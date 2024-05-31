function dcdt = Smad_model_trimer_screen_noise(t,c,k,PPase,Rin)

S2c = c(1);
S2n = c(2);
pS2c = c(3);
pS2n = c(4);
S4c = c(5);
S4n = c(6);
pS22c = c(7);
pS22n = c(8);
pS24c = c(9);
pS24n = c(10);
pS224c = c(11);
pS224n = c(12); %nuc. trimer

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

Ract = Rin(round(t,0)+1);

dS2c = -kphos*S2c*Ract -kin*S2c + kex*S2n;
dS2n = kdephos*pS2n*PPase + kin*S2c -kex*S2n;
dpS2c = kphos*S2c*Ract - 2*kon*pS2c^2 +  2*koff*pS22c - kon*pS2c*S4c +  koff*pS24c - kon*pS2c*pS24c + koff*pS224c - kin*pS2c + kex*pS2n;
dpS2n = -kdephos*pS2n*PPase - 2*kon*pS2n^2 + 2*koff*pS22n - kon*pS2n*S4n + koff*pS24n - kon*pS2n*pS24n + koff*pS224n + kin*pS2c - kex*pS2n;
dS4c = -kon*pS2c*S4c + koff*pS24c -kon*pS22c*S4c + koff*pS224c -kin*S4c + kin*S4n;
dS4n = -kon*pS2n*S4n + koff*pS24n - kon*pS22n*S4n + koff*pS224n + kin*S4c - kin*S4n;
dpS22c = kon*pS2c^2 - koff*pS22c - kon*pS22c*S4c + koff*pS224c - kin*pS22c +kex*pS22n;
dpS22n = kon*pS2n^2- koff*pS22n -kon*pS22n*S4n + koff*pS224n + kin*pS22c -kex*pS22n;
dpS24c = kon*pS2c*S4c - koff*pS24c - kon*pS24c*pS2c + koff*pS224c - kin*pS24c +kex*pS24n;
dpS24n = kon*pS2n*S4n -koff*pS24n -kon*pS24n*pS2n + koff*pS224n + kin*pS24c -kex*pS24n;
dpS224c = kon*pS2c*pS24c - koff*pS224c + kon*pS22c*S4c - koff*pS224c - CIF*kin*pS224c +kex*pS224n;
dpS224n = kon*pS2n*pS24n - koff*pS224n + kon*pS22n*S4n - koff*pS224n + CIF*kin*pS224c -kex*pS224n;

dcdt = [dS2c;dS2n;dpS2c;dpS2n;dS4c;dS4n;dpS22c;dpS22n;dpS24c;dpS24n;dpS224c;dpS224n];

end

