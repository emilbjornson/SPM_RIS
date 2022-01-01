function vn = refcoefficient(omega,Cn)
% Compute the reflection coefficient using the circuit from (3) in
% "Intelligent Reflecting Surface: Practical Phase Shift Model and
% Beamforming Optimization" by Samith Abeywickrama, Rui Zhang, Chau Yuen.

L1 = 2.5e-9; %Bottom layer inductance
L2 = 0.7e-9; %Top layer inductance
Rn = 1; %Effective resistance
Z0 = 377; %Impedance of free space

%Compute the impedance of the surface
Zn = 1i*omega*L1*(1i*omega*L2+1./(1i*omega*Cn)+Rn)./(1i*omega*L1+(1i*omega*L2+1./(1i*omega*Cn)+Rn));

%Compute reflection coefficient
vn = (Zn-Z0)./(Zn+Z0);

end