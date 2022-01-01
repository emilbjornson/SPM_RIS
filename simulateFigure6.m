%%This Matlab script generates parts of Figure 6 in the article:
%
%Emil Björnson, Henk Wymeersch, Bho Matthiesen, Petar Popovski, Luca
%Sanguinetti, and Elisabeth de Carvalho “Reconfigurable Intelligent
%Surfaces: A Signal Processing Perspective With Wireless Applications,”
%IEEE Signal Processing Magazine, To appear in March 2021.
%
%Download article: https://arxiv.org/pdf/2102.00742.pdf
%
%This is version 1.0 (Last edited: 2022-01-01)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%paper as described above.


close all;
clear;

%Figure 6(a) is generated with "true". Change to "false" to get 6(b).
A_LOS = true;


%% Define the array geometry

%Carrier frequency
fc = 3e9;

%Speed of light
c = 3e8;

%Wavelength
lambda = c/fc;

N_H = 20; %Number of elements per row in the RIS
N_V = 20; %Number of elements per column in the RIS
d_H = 0.25*lambda; %Horizontal element spacing
d_V = 0.25*lambda; %Vertical element spacing

N = N_H*N_V; %Total number of elements
U = zeros(3,N); %Matrix containing the position of the elements

i = @(m) mod(m-1,N_H); %Horizontal index
j = @(m) floor((m-1)/N_H); %Vertical index
for m = 1:N
    U(:,m) = [0; i(m)*d_H; j(m)*d_V]; %Position of the mth RIS element
end


%% Define the propagation environment

%Assume a LOS path between the RIS and UE
B_LOS = true;

%Location of transmitting AP (in meters)
APlocation = [40 -200 0];

%Location of receiving UE (in meters)
UElocation = [20 0 0]; 

distAP_RIS = norm(APlocation); %Distance from AP to the center of the RIS
distUE_RIS = norm(UElocation); %Distance from UE to the center of the RIS
distAP_UE = norm(APlocation - UElocation); %Distance from AP to UE

%Compute azimuth angles to the AP and UE as seen from the RIS
varphiAP_RIS = atan(APlocation(2)/APlocation(1));
varphiUE_RIS = atan(UElocation(2)/UElocation(1));

%Micro-cell pathloss models from ETSI TR 125 996 (where d is in meter)
pathlossLOS = @(d) db2pow(-30.18-26*log10(d));
pathlossNLOS = @(d) db2pow(-34.53-38*log10(d));


%Number of paths (including LOS paths)
La = 101;
Lb = 51;
Ld = 100;



%Subcarrier spacing
subSpacing = 15e3;

%Range of number of subcarriers
subCarriers = [25 50:50:1000];

%Compute the corresponding bandwidths
bandwidths = subSpacing*subCarriers;


%Number of Monte Carlo setups with random multipath components
realizations = 100;

%Prepare to compute the rates at the different methods
rateScatter = zeros(length(bandwidths),max(subCarriers),realizations);
rateNone = zeros(length(bandwidths),max(subCarriers),realizations);
rateLOSopt = zeros(length(bandwidths),max(subCarriers),realizations);
rateUpper = zeros(length(bandwidths),max(subCarriers),realizations);



%% Go through random realizations of the multipath components
for r = 1:realizations
    
    disp(['Realization ' num2str(r) ' out of ' num2str(realizations)]);
    
    
    %% Compute random propagation environment
    
    %Compute azimuth angles of arrival/departure of paths to/from the RIS. The
    %first path is LOS while other paths are distributed as in ETSI TR 125 996.
    deviationAzimuthAngleInterval = 40*pi/180; %Uniform distribution -40 to +40 degrees
    La_varphi = varphiAP_RIS + [0; deviationAzimuthAngleInterval*(2*rand(La-1,1)-1)];
    Lb_varphi = varphiUE_RIS + [0; deviationAzimuthAngleInterval*(2*rand(Lb-1,1)-1)];
    
    %Compute elevation angles of arrival/departure of paths to/from the RIS.
    %The first path is LOS while other paths are randomly distributed
    deviationElevationAngleInterval = 10*pi/180; %Uniform distribution -10 to +10 degrees
    La_theta = 0 + [0; deviationElevationAngleInterval*(2*rand(La-1,1)-1)];
    Lb_theta = 0 + [0; deviationElevationAngleInterval*(2*rand(Lb-1,1)-1)];
    
    %Propagation delays are computed based on the length of the LOS path.
    %Scattered paths have a delay that is uniformly distributed between the
    %LOS delay and twice this value.
    La_delay = distAP_RIS*(1 + [0; rand(La-1,1)])/c;
    Lb_delay = distUE_RIS*(1 + [0; rand(Lb-1,1)])/c;
    Ld_delay = distAP_UE*(1 + rand(Ld,1))/c;
    
    %The propagation loss of the different paths are computed so that the
    %LOS has its share according to the Rice factors, while the remaining
    %power is distributed among the other terms as in ETSI TR 125 996.
    La_powfactor = 10.^(-La_delay(2:end)+0.2*randn(La-1,1));
    Lb_powfactor = 10.^(-Lb_delay(2:end)+0.2*randn(Lb-1,1));
    Ld_powfactor = 10.^(-Ld_delay+0.2*randn(Ld,1));
    
    %The propagation losses from ETSI TR 125 996 are for isotropic
    %antennas. This factor compensates for the smaller size of the RIS
    %elements
    lossComparedToIsotropic = d_H*d_V/(lambda^2/(4*pi));
    
    %Compute propagation losses for channel to the RIS
    if A_LOS == true
        
        %Compute Rice factors according to ETSI TR 125 996
        RicefactorAP_RIS = db2pow(13-0.03*distAP_RIS);
        La_pathloss = lossComparedToIsotropic*pathlossLOS(distAP_RIS)*[RicefactorAP_RIS/(1+RicefactorAP_RIS); (1/(1+RicefactorAP_RIS))*La_powfactor/sum(La_powfactor)];
        
    else
        
        %Remove LOS path and sort the paths so the first one is strongest
        La = La-1;
        [La_pathloss,index] = sort(lossComparedToIsotropic*pathlossNLOS(distAP_RIS)*La_powfactor/sum(La_powfactor),'descend');
        La_varphi = La_varphi(index+1);
        La_theta = La_theta(index+1);
        La_delay = La_delay(index+1);
        
    end
    
    
    %Compute propagation losses for channel from the RIS
    if B_LOS == true
        
        %Compute Rice factors according to ETSI TR 125 996
        RicefactorUE_RIS = db2pow(13-0.03*distUE_RIS);
        Lb_pathloss = lossComparedToIsotropic*pathlossLOS(distUE_RIS)*[RicefactorUE_RIS/(1+RicefactorUE_RIS); (1/(1+RicefactorUE_RIS))*Lb_powfactor/sum(Lb_powfactor)];
        
    else
        
        %Remove LOS path and sort the paths so the first one is strongest
        Lb = Lb-1;
        [Lb_pathloss,index] = sort(lossComparedToIsotropic*pathlossNLOS(distUE_RIS)*Lb_powfactor/sum(Lb_powfactor),'descend');
        Lb_varphi = Lb_varphi(index+1);
        Lb_theta = Lb_theta(index+1);
        Lb_delay = Lb_delay(index+1);
        
    end
    
    
    %Compute propagation losses for the direct paths
    Ld_pathloss = pathlossNLOS(distAP_UE)*Ld_powfactor/sum(Ld_powfactor);
    
    
    
    %% Go through the range of considered bandwidths and compute discrete-time channels
    for b = 1:length(bandwidths)
        
        %Extract the bandwidth
        B = bandwidths(b);
        
        %Compute discrete-time impulse responses
        
        %Take the first sample when the first direct path arrives
        eta = min(Ld_delay);
        
        %Approximate the length of the impulse response based on the delay
        %spread plus 11 to compensate for the length of the pulse shape
        M = round(B*(2*(distAP_RIS+distUE_RIS)-distAP_UE)/c)+11;
        
        %Compute the number of subcarriers
        K = floor(B/subSpacing);

        %Compute V and h_d based on the formulas in the paper
        V = zeros(N,K);
        V_LOS = zeros(N,K);
        hd = zeros(K,1);
        
        %Go through all channel taps
        for k = 1:M
            
            %Go through all paths to and from the RIS
            for l1 = 1:La
                for l2 = 1:Lb
                    
                    %Compute elements of V based on Eq. (22)
                    V(:,k) = V(:,k) + sqrt(La_pathloss(l1)*Lb_pathloss(l2))*exp(-1i*2*pi*fc*(La_delay(l1)+Lb_delay(l2)))* functionSpatialSignature3DLoS(U,La_varphi(l1),La_theta(l1),lambda).* functionSpatialSignature3DLoS(U,Lb_varphi(l2),Lb_theta(l2),lambda)*sinc(k-1+B*(eta-La_delay(l1)-Lb_delay(l2)));
                    
                    %Extract only the strongest/LOS path, if it is used for RIS configuration
                    if (l1 == 1) && (l2 == 1)
                        V_LOS(:,k) = V(:,k);
                    end
                    
                end
            end
            
            %Go through all the direct paths
            for l = 1:Ld
                hd(k,1) = hd(k,1) + sqrt(Ld_pathloss(l))*exp(-1i*2*pi*fc*Ld_delay(l))*sinc(k-1+B*(eta-Ld_delay(l)));
            end
            
        end
        
        
        
        %% Compute achievable rates
        
        %Transmit power
        P = 1*B/1e6; %1 W per MHz
        
        %Noise power spectral density (including 10 dB noise figure)
        N0 = db2pow(-174+10)/1000; %W per Hz
        
        %Compute a K x K DFT matrix
        F = fft(eye(K));
        
        %Heuristic optimization of the RIS based on LOS path and the strongest tap
        [~,ind] = max(abs(hd)+sum(abs(V_LOS),1)');
        RIS_LOS_config = exp(1i*(angle(hd(ind))-angle(V_LOS(:,ind))));
        
        %Baseline with random configuration
        RIS_none_config = zeros(N,1);
        
        %Baseline with no delays configuration
        RIS_scattering = ones(N,1);
        

        %Compute the SNR per subcarrier that obtained with the different
        %RIS configurations if equal power allocation is used
        SNRsScatter = zeros(K,1);
        SNRsNone = zeros(K,1);
        SNRsLOSopt = zeros(K,1);
        SNRsUpper = zeros(K,1);
        
        for k = 1:K
            
            SNRsScatter(k) = P * abs(F(k,:)*(hd+V.'*RIS_scattering))^2/(B*N0);
            SNRsNone(k) = P * abs(F(k,:)*(hd+V.'*RIS_none_config))^2/(B*N0);
            SNRsLOSopt(k) = P * abs(F(k,:)*(hd+V.'*RIS_LOS_config))^2/(B*N0);
            SNRsUpper(k) = P * (abs(F(k,:)*hd) + sum(abs(F(k,:)*V.')))^2/(B*N0);
            
        end
        
        %Compute the rates that are obtained with the different RIS
        %configurations if optimal waterfilling power allocation is used
        powerAllocation = functionWaterfilling(K,1./SNRsScatter);
        rateScatter(b,1:K,r) =  B/(K+M-1)*log2(1+(SNRsScatter.*powerAllocation)');
        
        powerAllocation = functionWaterfilling(K,1./SNRsNone);
        rateNone(b,1:K,r) =  B/(K+M-1)*log2(1+SNRsNone.*powerAllocation);
        
        powerAllocation = functionWaterfilling(K,1./SNRsLOSopt);
        rateLOSopt(b,1:K,r) =  B/(K+M-1)*log2(1+SNRsLOSopt.*powerAllocation);
        
        powerAllocation = functionWaterfilling(K,1./SNRsUpper);
        rateUpper(b,1:K,r) =  B/(K+M-1)*log2(1+SNRsUpper.*powerAllocation);
        
    end
    
    
end


set(groot,'defaultAxesTickLabelInterpreter','latex');


%% Plot the simulation results
figure;
hold on; box on; grid on;
plot(bandwidths/1e6,mean(sum(rateUpper,2),3)/1e6,'r-','LineWidth',2);
plot(bandwidths/1e6,mean(sum(rateLOSopt,2),3)/1e6,'k--','LineWidth',2);
plot(bandwidths/1e6,mean(sum(rateScatter,2),3)/1e6,'k:','LineWidth',2);
plot(bandwidths/1e6,mean(sum(rateNone,2),3)/1e6,'b-.','LineWidth',2);
set(gca,'fontsize',16);
xlabel('Bandwidth [MHz]','Interpreter','Latex');
ylabel('Rate [Mbps]','Interpreter','Latex');
legend({'Upper bound','Heuristic: STM','No RIS: Metal sheet','No RIS: Absorption'},'Interpreter','latex','Location','NorthWest');
