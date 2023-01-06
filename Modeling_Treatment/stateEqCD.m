%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State equation for CD
%
% Author: Tuğba Akman Date: Jan 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dy = stateEqCD(t,y,uCD,Tu)

global Eq

dy = zeros(4,1);

uCD = interp1(Tu,uCD,t); 

SCD = y(1);
RCD = y(2);
ECD = y(3);
FCD = y(4);

dy(1) = Eq.k1*SCD*(1- ((SCD+Eq.m*RCD)/Eq.K))*((ECD^Eq.n)/(Eq.a1+(ECD^Eq.n)))- Eq.betaCD*((Eq.a3^Eq.l)/(Eq.a3^Eq.l + ECD^Eq.l))*SCD- ((Eq.a2^Eq.l)/(Eq.a2^Eq.l + ECD^Eq.l))*SCD;
dy(2) = Eq.k3CD*RCD*(1- ((SCD+Eq.m*RCD)/Eq.K)) + Eq.betaCD*((Eq.a3^Eq.l)/(Eq.a3^Eq.l + ECD^Eq.l))*SCD;
dy(3) = (1-uCD)*Eq.pCD*Eq.r*FCD - Eq.mu*ECD ;
dy(4) = Eq.k2CD*FCD*(1- (FCD*Eq.m2CD))- Eq.alpha*(SCD+RCD)*FCD;


end