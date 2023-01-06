%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State equation for HFD
%
% Author: Tuğba Akman Date: Jan 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dy = stateEqHFD(t,y,uHFD,Tu)

global Eq

dy = zeros(4,1);

uHFD = interp1(Tu,uHFD,t); 

SHFD = y(1);
RHFD = y(2);
EHFD = y(3);
FHFD = y(4);

dy(1) = Eq.k1*SHFD*(1- ((SHFD+Eq.m*RHFD)/Eq.K))*((EHFD^Eq.n)/(Eq.a1+(EHFD^Eq.n)))- Eq.betaHFD*((Eq.a3^Eq.l)/(Eq.a3^Eq.l + EHFD^Eq.l))*SHFD-  ((Eq.a2^Eq.l)/(Eq.a2^Eq.l + EHFD^Eq.l))*SHFD;
dy(2) = Eq.k3HFD*RHFD*(1- ((SHFD+Eq.m*RHFD)/Eq.K)) + Eq.betaHFD*((Eq.a3^Eq.l)/(Eq.a3^Eq.l + EHFD^Eq.l))*SHFD;
dy(3) = (1-uHFD)*Eq.pHFD*Eq.r*FHFD - Eq.mu*EHFD;
dy(4) = Eq.k2HFD*FHFD*(1- (FHFD*Eq.m2HFD))- Eq.alpha*(SHFD+RHFD)*FHFD;


end