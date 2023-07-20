function [anSOP] = Correlated_SOP_UserN(SNR,aSN,aSE,muSN_dB,muSE_dB,sigmaSN_dB,sigmaSE_dB,eta,rS,bN,pS,rho)

%This function calculates the analytic result of the SOP of the near user
%under external eavesdropping using the Gauss-Chebyshev quadrature

thetaN = 2^rS;

temp1 = 0;
K = 100000;






muX_dB = 2*muSN_dB; muX= (0.1*(log(10))*muX_dB);
muZ_dB = 2*muSE_dB; muZ= (0.1*(log(10))*muZ_dB);
sigX_dB = 2*sigmaSN_dB; sigX= (0.1*(log(10))*sigX_dB);
sigZ_dB = 2*sigmaSE_dB; sigZ= (0.1*(log(10))*sigZ_dB);

phi = (pS*aSN^2)/SNR;
varrho = (pS*aSE^2)/SNR;

for k = 1: K

    vartheta = pi/K;
    vk = cos(((2*k - 1)*pi)/(2*K));
    nu = (pi/4)*(vk + 1);
    psi = tan(nu);
    
    temp1 = temp1 + ((eta*pi*vartheta*sqrt(1-vk^2)*sec(nu)^2)/(8*psi*sqrt(2*pi)*sigZ_dB))* ...
           exp(-((eta*log(psi) - muZ_dB)^2)/((2*sigZ_dB^2)))* ...
          ((1 + erf((((eta*log(thetaN*bN*varrho*psi + thetaN - 1) - eta*log(bN*phi) - muX_dB)/(sqrt(2)*sigX_dB)) - rho*((eta*log(psi) - muZ_dB )/(sqrt(2)*sigX_dB)))/(sqrt(1 - rho^2)))));
end
   
anSOP = temp1;
end

