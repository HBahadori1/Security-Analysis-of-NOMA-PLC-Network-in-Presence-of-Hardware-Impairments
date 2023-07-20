function [anSOP] = Correlated_SOP_UserF(SNR,aSF,aSE,muSF_dB,muSE_dB,sigmaSF_dB,sigmaSE_dB,eta,rS,bN,bF,pS,rho)

%This function calculates the analytic result of the SOP of the far user
%under external eavesdropping using the Gauss-Chebyshev quadrature

thetaF = 2^rS;

temp1 = 0;
K = 100000;



muY_dB = 2*muSF_dB ; muY= (0.1*(log(10))*muY_dB);
muZ_dB = 2*muSE_dB ; muZ= (0.1*(log(10))*muZ_dB);
sigY_dB = 2*sigmaSF_dB; sigY= (0.1*(log(10))*sigY_dB);
sigZ_dB = 2*sigmaSE_dB; sigZ= (0.1*(log(10))*sigZ_dB);

beta = (pS*aSF^2)/SNR;
varrho = (pS*aSE^2)/SNR;

for k = 1: K
    
    UL = (bF - bN*(thetaF-1))/(bN*varrho*(bF+bN)*(thetaF -1));
    w = pi/K;
    vk = cos(((2*k - 1)*pi)/(2*K));
    psi = (UL/2)*(vk + 1);
    
    lambda = (thetaF-1)*bN*varrho + thetaF*bF*varrho;
    varpi = bN*varrho*(bF+bN)*(thetaF-1);
    
    
    PartA = (lambda*psi +(thetaF -1))/(beta*(bF-bN*(thetaF-1) - varpi*psi));
    
    
    temp1 = temp1 + ((eta*UL*w*sqrt(1-vk^2))/(2*psi*sqrt(2*pi)*sigZ_dB))* ...
           exp(-((eta*log(psi) - muZ_dB)^2)/((2*sigZ_dB^2)))* ...
           (1-((1/2)*(1 + erf((((eta*log(PartA) - muY_dB)/(sqrt(2)*sigY_dB)) - rho*((eta*log(psi) - muZ_dB)/(sqrt(2)*sigZ_dB)))/(sqrt(1 - rho^2))))));
end


   
anSOP = 1 - temp1;
end
