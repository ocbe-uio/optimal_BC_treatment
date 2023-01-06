%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CoState equation for HFD
%
% Author: Tuğba Akman Date: Jan 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dp = costateEqHFD(t,p,uHFD,Tu,x5,x6,x7,x8,xt)

global Eq

dp = zeros(4,1);

x5 = interp1(xt,x5,t); 
x6 = interp1(xt,x6,t); 
x7 = interp1(xt,x7,t); 
x8 = interp1(xt,x8,t); 

uHFD = interp1(Tu,uHFD,t); 

dp(1) = -Eq.weight4 - p(1).*((Eq.k1*(x7.^Eq.n)/(Eq.a1 + x7.^Eq.n))*(1-((2*x5+Eq.m*x6)/Eq.K)) - Eq.betaHFD*((Eq.a3^Eq.l)./(Eq.a3^Eq.l + x7.^Eq.l)) - ((Eq.a2^Eq.l)./(Eq.a2^Eq.l + x7.^Eq.l)))...
      - p(2).*(-(Eq.k3HFD*(x6/Eq.K)) + Eq.betaHFD.*((Eq.a3^Eq.l)./(Eq.a3^Eq.l + x7.^Eq.l)))... 
      + p(4).*(Eq.alpha*x8) ;

dp(2) =  -Eq.weight5- p(1).*((Eq.k1*(x7.^Eq.n)/(Eq.a1 + x7.^Eq.n)).*x5*(-Eq.m/Eq.K))...
      - p(2).*(Eq.k3HFD*(1-((x5+2*Eq.m*x6)/Eq.K)))...
      + p(4).*(Eq.alpha*x8);
  
dp(3) = -p(1).*((((Eq.k1*Eq.a1*Eq.n*(x7.^(Eq.n-1)))./(Eq.a1+x7.^Eq.n).^2)).*x5.*(1-((x5+Eq.m*x6)/Eq.K)) + Eq.betaHFD*((Eq.a3^Eq.l)./((Eq.a3^Eq.l + x7^Eq.l).^2)).*(Eq.l).*(x7.^(Eq.l-1)).*x5...
      +((Eq.a2^Eq.l)./((Eq.a2^Eq.l + x7^Eq.l).^2)).*(Eq.l).*(x7.^(Eq.l-1)).*x5) ...
      + p(2).*Eq.betaHFD.*((Eq.a3^Eq.l)./((Eq.a3^Eq.l + x7^Eq.l).^2)).*(Eq.l).*(x7.^(Eq.l-1)).*x5...  
      + p(3).*Eq.mu;
  
dp(4) = -p(3).*(Eq.pHFD*Eq.r*(1-uHFD)) - p(4).*(Eq.k2HFD - 2*Eq.k2HFD*Eq.m2HFD*x8 - Eq.alpha*(x5+x6));


end