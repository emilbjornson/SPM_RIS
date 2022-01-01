%%This Matlab script generates parts of Figure 4 in the article:
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

%Set propagation losses per element
alpha = db2pow(-80); %To RIS
beta = db2pow(-60); %From RIS
rho1 = db2pow(-80); %Option 1 for uncontrollable path
rho2 = db2pow(-110); %Option 2 for uncontrollable path
gamma = 1; %Loss in RIS

%Bandwidth (Hz)
B = 1e6;

%Number of RIS elements
N = 0:400;

%Set transmit SNR (dependable on bandwidth)
PBN0 = db2pow(100);

%Compute SNRs using (31) with ideal RIS configuration
SNR1 = PBN0*(sqrt(rho1)+N*sqrt(alpha*beta*gamma)).^2;
SNR2 = PBN0*(sqrt(rho2)+N*sqrt(alpha*beta*gamma)).^2;


%Number of realizations in Monte Carlo simulation
numberOfRealizations = 10000;

%The coordinate system is selected so that the uncontrollable path has zero
%phase, while all other paths have uniformly distributed phases
phaseShifts = rand(max(N),numberOfRealizations)*2*pi;

%The phase of each path is rotated to be between -pi/4 and pi/4, using an
%RIS with 4 configurations
rotateShifts = mod(phaseShifts,pi/2)-pi/4;

%Compute the actual SNRs with an RIS with 4 configurations
SNR1b = PBN0*[rho1; mean(abs(sqrt(rho1)+cumsum(exp(1i*rotateShifts),1)*sqrt(alpha*beta*gamma)).^2,2)];
SNR2b = PBN0*[rho2; mean(abs(sqrt(rho2)+cumsum(exp(1i*rotateShifts),1)*sqrt(alpha*beta*gamma)).^2,2)];


set(groot,'defaultAxesTickLabelInterpreter','latex');


%% Plot the simulation results
figure;
hold on; box on; grid on;
plot(N,B*log2(1+SNR1)/1e6,'r-','LineWidth',2)
plot(N,B*log2(1+SNR1b)/1e6,'k--','LineWidth',2)
plot(N,B*log2(1+SNR2)/1e6,'b-.','LineWidth',2)
plot(N,B*log2(1+SNR2b)/1e6,'k:','LineWidth',2)
set(gca,'fontsize',16);
xlabel('Number of RIS elements ($N$)','Interpreter','latex');
ylabel('Rate [Mbps]','Interpreter','latex');
legend({'$\rho=-80$ dB (ideal config)','$\rho=-80$ dB (4 configs)','$\rho=-110$ dB (ideal config)','$\rho=-110$ dB (4 configs)'},'Interpreter','latex','Location','NorthWest','Location','Best');

