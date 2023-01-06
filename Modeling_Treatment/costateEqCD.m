%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CoState equation for CD
%
% Author: Tuğba Akman Date: Jan 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dp = costateEqCD(t,p,uCD,Tu,x1,x2,x3,x4,xt)

global Eq

dp = zeros(4,1);

x1 = interp1(xt,x1,t); 
x2 = interp1(xt,x2,t);
x3 = interp1(xt,x3,t); 
x4 = interp1(xt,x4,t); 
 

uCD = interp1(Tu,uCD,t); 

dp(1) = -Eq.weight1 - p(1).*((Eq.k1*(x3.^Eq.n)/(Eq.a1 + x3.^Eq.n))*(1-((2*x1+Eq.m*x2)/Eq.K)) - Eq.betaCD*((Eq.a3^Eq.l)./(Eq.a3^Eq.l + x3.^Eq.l)) - ((Eq.a2^Eq.l)./(Eq.a2^Eq.l + x3.^Eq.l)))...
      - p(2).*(-(Eq.k3CD*(x2/Eq.K)) + Eq.betaCD*((Eq.a3^Eq.l)./(Eq.a3^Eq.l + x3.^Eq.l)))... 
      + p(4).*(Eq.alpha*x4) ;

dp(2) =  -Eq.weight2- p(1).*((Eq.k1*(x3.^Eq.n)/(Eq.a1 + x3.^Eq.n)).*x1*(-Eq.m/Eq.K))...
      - p(2).*(Eq.k3CD*(1-((x1+2*Eq.m*x2)/Eq.K)))...
      + p(4).*(Eq.alpha*x4);

dp(3) = -p(1).*((((Eq.k1*Eq.a1*Eq.n*(x3.^(Eq.n-1)))./(Eq.a1+x3.^Eq.n).^2)).*x1.*(1-((x1+Eq.m*x2)/Eq.K)) + Eq.betaCD*((Eq.a3^Eq.l)./((Eq.a3^Eq.l + x3^Eq.l).^2)).*(Eq.l).*(x3.^(Eq.l-1)).*x1...
      +((Eq.a2^Eq.l)./((Eq.a2^Eq.l + x3^Eq.l).^2)).*(Eq.l).*(x3.^(Eq.l-1)).*x1) ...
      + p(2).*Eq.betaCD.*((Eq.a3^Eq.l)./((Eq.a3^Eq.l + x3^Eq.l).^2)).*(Eq.l).*(x3.^(Eq.l-1)).*x1...  
      + p(3).*Eq.mu;
  
dp(4) = -p(3).*(Eq.pCD*Eq.r*(1-uCD)) - p(4).*(Eq.k2CD - 2*Eq.k2CD*Eq.m2CD*x4 - Eq.alpha*(x1+x2));
  


end