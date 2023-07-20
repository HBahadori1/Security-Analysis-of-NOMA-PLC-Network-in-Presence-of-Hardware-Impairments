%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%Physical Layer Security in NOMA PLC%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%Hamidreza Bahadori December 2022%%%%%%%%%%%%%%%%%%%%%%

%This matlab code is used to obtain the SOP of the near and far users under
%internal and external eavesdropping, respectively.
%Correlation factor, rho varies between 0 and 1
%rho = 0 denotes independent channels
%rho = 1 denotes fully correlated channels

% -------------------------------------------------------------------------
% -------------------------- INITIAL VARIABLES ----------------------------
% -------------------------------------------------------------------------
clear all
close all
clc

SNR_dB = -10:3:30; %SBNR of the system signal of background noise
SNR = 10.^( (-SNR_dB)/10);% conversion of SNR.
SINR_dB = -10; %%SINR of the system signal of impulsive noise
SINR = 10.^( (-SINR_dB)/10); 
pS = 1; %the transmit power at source in watts
dSN = 80; %distance of S-->N link.
dSF = 120; %distance of S-->F link.
dSE = 150; %distance of S-->E link.

b0 = 9.4*1e-3; %attenuation constant
b1 = 4.2*1e-7; %attenuation constant
f = 30; %frequency
k = 0.7; %attenuation constant
eta=10/log(10);
p = 0.001; %probability of occurrence of impulsive noise


aSN = exp(-(b0+b1*(f^k))*dSN); %S-->N link attenuation.
aSF = exp(-(b0+b1*(f^k))*dSF); %S-->F link attenuation.
aSE = exp(-(b0+b1*(f^k))*dSE); %S-->E link attenuation.

muSN_dB = 3; muSN= (0.1*(log(10))*muSN_dB); %channel mean of S-->N link.
muSF_dB = 3; muSF= (0.1*(log(10))*muSF_dB); %channel mean of S-->F link.
muSE_dB = 3; muSE= (0.1*(log(10))*muSE_dB); %channel mean of S-->E link.
sigmaSN_dB = 1; sigmaSN= (0.1*(log(10))*sigmaSN_dB); %standard deviation of S-->N link.
sigmaSF_dB = 1; sigmaSF= (0.1*(log(10))*sigmaSF_dB); %standard deviation of S-->F link.
sigmaSE_dB = 1; sigmaSE= (0.1*(log(10))*sigmaSE_dB); %standard deviation of S-->E link.

%Correlation Paramters
Mu = [muSN muSF muSE];
Sigma = [sigmaSN sigmaSF sigmaSE];
rho = 0.2; %correlation factor varies between 0 and 1
corrM = [1 rho rho; 
         rho 1 rho; 
         rho rho 1];
Simulations = 1;

%Power allocation
bN = 0.1; %power allocation of S-->N link.
bF = 0.9; %power allocation of S-->F link.

%Secrecy target rates
rN = 0.1; % near user target secrecy rate
rF = 0.1; % far user target secrecy rate



numchannel = 100000; %number of channel realizations
nSOP_an = zeros(length(p), length(SNR)); %store analytic results
fSOP_an = zeros(length(p), length(SNR)); %store analytic results


out_event_N_BG =  zeros(length(p), length(SNR)); %store sim results
out_event_N_IN =  zeros(length(p), length(SNR)); %store sim results

out_event_F_BG =  zeros(length(p), length(SNR)); %store sim results
out_event_F_IN =  zeros(length(p), length(SNR)); %store sim results


for pcount = 1: length(p)

for idx_SNR = 1: length(SNR)
    
    for idx_ch = 1: numchannel% channel begin
        
        [pcount idx_SNR idx_ch];
        
         %Generate correlated channels 
         y = MvLogNRand( Mu, Sigma, Simulations , corrM );
         hSN = y(1); %source to near user channel
         hSF = y(2); %source to far user channel 
         hSE = y(3); %source to eve channel
        
        
        % -----------------------------------------------------------------
        % ------------------EXTERNAL EAVESDROPPER NEAR USER----------------
        % -----------------------------------------------------------------
        yNxn_BG = ((bN*pS*aSN^2*hSN^2)/(SNR(idx_SNR )));  %SNR of near user with BG noise
        yNxn_IN = ((bN*pS*aSN^2*hSN^2)/(SNR(idx_SNR ) + SINR)); %SNR of near user with BG and IN noise
        
        yExn_BG = ((bN*pS*aSE^2*hSE^2)/(SNR(idx_SNR ))) ; %SNR of eve with BG noise
        yExn_IN = ((bN*pS*aSE^2*hSE^2)/(SNR(idx_SNR ) + SINR)); %SNR of eve with BG and IN noise
        
        
        cNxn_BG = log2(1 + yNxn_BG); %cap of near user with BG noise
        cExn_BG = log2(1 + yExn_BG); %cap of eve with BG noise
        cS_xn_BG = max(0, cNxn_BG -cExn_BG); %secrecy cap of near user with BG noise
        
        cNxn_IN = log2(1 + yNxn_IN); %cap of near user with BG and IN noise
        cExn_IN = log2(1 + yExn_IN); %cap of eve with BG and IN noise 
        cS_xn_IN = max(0, cNxn_IN -cExn_IN); %secrecy cap of near user with BG and IN noise
        
        
        if cS_xn_BG < rN
            out_event_N_BG(pcount,idx_SNR) = out_event_N_BG(pcount,idx_SNR) + 1; %outage event for BG noise
        end
         
        if cS_xn_IN < rN
            out_event_N_IN(pcount,idx_SNR) = out_event_N_IN(pcount,idx_SNR) + 1; %outage event for BG and IN noise 
        end
        
        

        
        %-----------------------------------------------------------------
        %------------------EXTERNAL EAVESDROPPER-FAR USER-----------------
        %-----------------------------------------------------------------
        yFxf_BG = ((bF*pS*aSF^2*hSF^2)/(bN*pS*aSF^2*hSF^2 + SNR(idx_SNR ))); %SNR of far user with BG noise
        yFxf_IN = ((bF*pS*aSF^2*hSF^2)/(bN*pS*aSF^2*hSF^2 + SNR(idx_SNR) + SINR)); %SNR of far user with BG and IN noise
        
        yExf_BG = ((bF*pS*aSE^2*hSE^2)/(bN*pS*aSE^2*hSE^2 + SNR(idx_SNR ))); %SNR of eve with BG noise
        yExf_IN = ((bF*pS*aSE^2*hSE^2)/(bN*pS*aSE^2*hSE^2 + SNR(idx_SNR) + SINR));  %SNR of eve with BG and IN noise
        
        cFxf_BG = log2(1 + yFxf_BG); %cap of far user with BG noise
        cExf_BG = log2(1 + yExf_BG); %cap of eve with BG noise
        cS_xf_BG = max(0, cFxf_BG - cExf_BG); %secrecy cap of far user with BG noise
        
        cFxf_IN = log2(1 + yFxf_IN); %cap of far user with BG and IN noise
        cExf_IN = log2(1 + yExf_IN); %cap of eve with BG and IN noise
        cS_xf_IN = max(0, cFxf_IN  -cExf_IN); %secrecy cap of far user with BG and IN noise
        
        if cS_xf_BG < rF
            out_event_F_BG(pcount,idx_SNR) = out_event_F_BG(pcount,idx_SNR) + 1;  %outage event for BG noise
        end
         
        if cS_xf_IN < rF
            out_event_F_IN(pcount,idx_SNR) = out_event_F_IN(pcount,idx_SNR) + 1;  %outage event for BG and IN noise
        end
        
        
          
 
    end
    
    
        % -------------------------------------------------------------------------
        % -----------------------------ANALYTIC RESULTS----------------------------
        % -------------------------------------------------------------------------
        %SOP of near user under external eavesdropping
        corr_userN_BG = Correlated_SOP_UserN(SNR(idx_SNR),aSN,aSE,muSN_dB,muSE_dB,sigmaSN_dB,sigmaSE_dB,eta,rN,bN,pS,rho);
        corr_userN_IN = Correlated_SOP_UserN(SNR(idx_SNR)+SINR,aSN,aSE,muSN_dB,muSE_dB,sigmaSN_dB,sigmaSE_dB,eta,rN,bN,pS,rho);
        nSOP_an(pcount,idx_SNR) = (1- p(pcount))*corr_userN_BG  + p(pcount)*corr_userN_IN;
        
 
        %SOP of far user under external eavesdropping
        corr_userF_BG = Correlated_SOP_UserF(SNR(idx_SNR),aSF,aSE,muSF_dB,muSE_dB,sigmaSF_dB,sigmaSE_dB,eta,rF,bN,bF,pS,rho);
        corr_userF_IN = Correlated_SOP_UserF(SNR(idx_SNR) +SINR,aSF,aSE,muSF_dB,muSE_dB,sigmaSF_dB,sigmaSE_dB,eta,rF,bN,bF,pS,rho);
        fSOP_an(pcount,idx_SNR) = (1- p(pcount))*corr_userF_BG + p(pcount)*corr_userF_IN;
        
end
end

nSOP = (1-p)*(out_event_N_BG/numchannel) + p*(out_event_N_IN/numchannel); %SOP of near user simulation
fSOP = (1-p)*(out_event_F_BG/numchannel) + p*(out_event_F_IN/numchannel); %SOP of far user simulation 
SOP_sys = (1 - (1 - nSOP).*(1 - fSOP)); % SOP of the system




figure(1)
semilogy(SNR_dB,nSOP_an,'k','LineWidth',1.5);
semilogy(SNR_dB,SOP_sys,'b--','LineWidth',1.5);
grid on
hold on
semilogy(SNR_dB,fSOP_an,'r-','LineWidth',1.5);
semilogy(SNR_dB,nSOP,'-ko');
semilogy(SNR_dB,fSOP,'-bo');
semilogy(SNR_dB,SOP_sys,'-ro');
%ylim([10^-2 1])
xlabel('SNR [dB]','Interpreter','Latex')
ylabel('SOP','Interpreter','Latex')
legend('System SOP','Near user','Far user','Simulation','Location','southwest')
 