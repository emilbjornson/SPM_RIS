%%This Matlab script generates parts of Figure 8 in the article:
%
%Emil Björnson, Henk Wymeersch, Bho Matthiesen, Petar Popovski, Luca
%Sanguinetti, and Elisabeth de Carvalho “Reconfigurable Intelligent
%Surfaces: A Signal Processing Perspective With Wireless Applications,”
%IEEE Signal Processing Magazine, To appear in March 2021.
%
%Download article: https://arxiv.org/pdf/2102.00742.pdf
%
%This is version 1.0 (Last edited: 2022-01-10)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%paper as described above.
%
% (c) 2021, Henk Wymeersch, henkw@chalmers.se

close all;
clear all;
warning off;
clc;

% 1. signal parameters
% ---------------------
signal.fc=28;                           % carrier in GHz
signal.c=0.3;                           % speed of light [m/ns]
signal.T=16;                            % number of OFDM symbols
signal.P_dBm=20;                        % transmit power dBm        
signal.P=10^(0.1*signal.P_dBm);         % transmit power mW
signal.N0dBmHz=-174;                    % noise PDS [dBm/Hz]
signal.N0=10^(0.1*signal.N0dBmHz)*1e9;  % noise PSD  [mW/GHz] (290 Kelvin * Boltzmann constant in W/Hz)
signal.BW=0.4;                          % Bandwidth GHz
signal.NFdB=0;                          % receiver noise figure [dB]
signal.Ns=3000;                         % active subcarriers
signal.Deltaf=signal.BW/signal.Ns;      % subcarrier spacing [GHz]
signal.Es=signal.P/(signal.Ns*signal.Deltaf);       % energy per subcarrier
signal.NF=10^(0.1*signal.NFdB);         % noise figure in linear scale
signal.sigma2=signal.NF*signal.N0/2;    % noise variance per subcarrier for real and imaginary component
signal.lambda=signal.c/signal.fc;       % wavelength
signal.data=ones(signal.T,signal.Ns);   % pilots (set to all ones)
signal.duration=signal.Ns*signal.T/signal.BW/1000000; % duration of the signal in milliseconds


% 2. default RIS, BS, UE, SP parameters
% -------------------------------------
RIS.M=64;                               % RIS elements (note: linear RIS)
RIS.Delta=signal.lambda/5;              % RIS element spacing
RIS.location=[-5;5];                    % location of RIS in XY plane
RIS.kappa=0;                            % RIS rotation around the Y axis 
RIS.size=RIS.Delta*RIS.M;               % total RIS size in meters
BS.location=[0; 0];                     % BS location in meters
UE.location=[-3;8];                     % user location in meters
SP.location=[3;4];                      % scatter point location in meters
SP.RCS=db2pow(0);                       % scatter point radar cross section


% 3. Offline Optimization of RIS placement with 1 or 3 RIS
% --------------------------------------------------------
x_grid=linspace(-7.5,7.5,50);
y_grid=linspace(-2,13,50);
plotPEBContour(signal,RIS,BS,[-5 0  5; 5  10 5; 0 +pi/2 +pi],x_grid,y_grid);    % PEB with 3 RIS
plotPEBContour(signal,RIS,BS,[0;10;+pi/2],x_grid,y_grid);                       % PEB with 1 RIS on wall facing the BS
plotPEBContour(signal,RIS,BS,[-5;5;0],x_grid,y_grid);                           % PEB with 1 RIS on left wall
plotPEBContour(signal,RIS,BS,[5;5;+pi],x_grid,y_grid);                          % PEB with 1 RIS on right wall

% 4. Online Optimization of RIS beams
% --------------------------------------------------------
onlinedesign(BS,UE,RIS,signal);    

disp('done!')