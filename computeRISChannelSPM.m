function [alpha,response,derivativeResponse,theta]=computeRISChannelSPM(location,RIS,signal)
% function [alpha,response,derivativeResponse,theta]=computeRISChannelSPM(location,RIS,signal)
% (c) 2021, Henk Wymeersch, henkw@chalmers.se
% This is version 1.0 (Last edited: 2022-01-10)
% 
% Operation: generates the channel gain, the response vector, the derivative of the
% response, and the angle of arrival
% inputs:
%   -location: a transmitter location
%   -RIS: RIS structure
%   -signal: signal structure
% 
% outputs
%   -alpha: complex channel gain
%   -response: Mx1 response vector
%   -derivativeResponse: Mx1 derivative of the response vector
%   -theta: angle of arrvial
    
    d = norm(location-RIS.location);
    R=[cos(RIS.kappa) -sin(RIS.kappa); sin(RIS.kappa) cos(RIS.kappa)];
    posR=R*(location-RIS.location);
    theta=atan2(posR(2),posR(1));
    ii=0:RIS.M-1;
    ii=ii-RIS.M/2+0.5; % center around zero
    a=exp(-1j*2*pi*RIS.Delta/signal.lambda*ii*sin(theta));
    a=a.';
    da=diag(-1j*2*pi*RIS.Delta/signal.lambda*ii*cos(theta))*a;    
    if (abs(theta)>pi/2)
        directivity=0;
    else        
        directivity=1;
    end                
    ar=(RIS.Delta);   % basic element size    
    alpha=exp(1j*rand(1)*2*pi)*directivity*ar/(sqrt(4*pi)*d);    % channel amplitude
    response=a;
    derivativeResponse=da;
