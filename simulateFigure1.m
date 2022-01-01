%%This Matlab script generates parts of Figure 1 in the article:
%
%Emil Björnson, Henk Wymeersch, Bho Matthiesen, Petar Popovski, Luca
%Sanguinetti, and Elisabeth de Carvalho “Reconfigurable Intelligent
%Surfaces: A Signal Processing Perspective With Wireless Applications,”
%IEEE Signal Processing Magazine, To appear in March 2022.
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


%Range of frequency values
fRange = linspace(2,4,100);

%Varying capacitance value
Cn = [0.79 0.88 0.96 2.2]*1e-12;

%Prepare to store results
vn = zeros(length(fRange),length(Cn));

for k = 1:length(fRange)

    %Compute angular frequency
    omega = 2*pi*fRange(k)*1e9;

    for m = 1:length(Cn)

        %Compute reflection coefficient using circuit from (3) in
        % "Intelligent Reflecting Surface: Practical Phase Shift Model and
        % Beamforming Optimization" by Samith Abeywickrama, Rui Zhang, Chau Yuen.
        vn(k,m) =  refcoefficient(omega,Cn(m));

    end

end


set(groot,'defaultAxesTickLabelInterpreter','latex');


%% Plot the simulation results
figure;
hold on; box on; grid on;
plot(fRange,angle(vn(:,1)),'k','LineWidth',2);
plot(fRange,angle(vn(:,2)),'r--','LineWidth',2);
plot(fRange,angle(vn(:,3)),'b-.','LineWidth',2);
plot(fRange,angle(vn(:,4)),'k:','LineWidth',2);
set(gca,'fontsize',16);
xlabel('Frequency ($f$) [GHz]','Interpreter','Latex');
ylabel('Phase response [rad]','Interpreter','Latex');
legend({[num2str(Cn(1)*1e12) ' pF'],[num2str(Cn(2)*1e12) ' pF'],[num2str(Cn(3)*1e12) ' pF'],[num2str(Cn(4)*1e12) ' pF']},'Interpreter','Latex','Location','NorthEast');
ylim([-pi pi]);
yticks(-pi:pi/2:pi);
yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})


figure;
hold on; box on; grid on;
plot(fRange,pow2db(abs(vn(:,1))),'k','LineWidth',2);
plot(fRange,pow2db(abs(vn(:,2))),'r--','LineWidth',2);
plot(fRange,pow2db(abs(vn(:,3))),'b-.','LineWidth',2);
plot(fRange,pow2db(abs(vn(:,4))),'k:','LineWidth',2);
set(gca,'fontsize',16);
xlabel('Frequency ($f$) [GHz]','Interpreter','Latex');
ylabel('Amplitude response [dB]','Interpreter','Latex');
legend({[num2str(Cn(1)*1e12) ' pF'],[num2str(Cn(2)*1e12) ' pF'],[num2str(Cn(3)*1e12) ' pF'],[num2str(Cn(4)*1e12) ' pF']},'Interpreter','Latex','Location','SouthEast');
ylim([-3 0]);
