function [PEB,J] = computePEBSPM(BS,UE,RIS,signal,beams_ref)
% function [PEB,J] = computePEBSPM(BS,UE,RIS,signal,beams_ref)
% (c) 2020, Henk Wymeersch, henkw@chalmers.se
% This is version 1.0 (Last edited: 2022-01-10)
% 
% Operation: computes the PEB and channel FIM for a system with LOS link
% and a path from 1 RIS
% inputs:
%   -BS: BS structure
%   -UE: UE structure
%   -RIS: RIS structure
%   -signal: signal structure
%   -beams_ref: a collection of beams
% 
% outputs:
%   -PEB: scalar position error bound in meters
%   -FIM: 7x7 FIM of the channel parameters, ordered as [tauLOS, tauNLOS, theta, R{alphaLOS}, I{alphaLOS}, R{alphaNLOS}, I{alphaNLOS}]

    % determine geometric parameters
    [alphaUE,aUE,daUE]=computeRISChannelSPM(UE.location,RIS,signal);
    [alphaBS,aBS,~]=computeRISChannelSPM(BS.location,RIS,signal);
    alphaLOS=signal.lambda/(4*pi*norm(UE.location-BS.location))*exp(1j*rand(1)*2*pi);   
    alphaNLOS=alphaUE*alphaBS*exp(1j*rand(1)*2*pi);
    tauLOS=norm(UE.location-BS.location)/signal.c;
    tauNLOS=norm(UE.location-RIS.location)/signal.c+norm(BS.location-RIS.location)/signal.c;
    ii=0:signal.Ns-1;
    ii=ii-signal.Ns/2+0.5; % cennter at zero
    dtauLOS= exp(-1j*2*pi*ii*signal.Deltaf*tauLOS);    
    dtauNLOS=exp(-1j*2*pi*ii*signal.Deltaf*tauNLOS);       
    ddtauLOS=(-1j*2*pi*ii*signal.Deltaf).*dtauLOS;    
    ddtauNLOS=(-1j*2*pi*ii*signal.Deltaf).*dtauNLOS;
    
    
    % compute FIM of channel parameters    
    Jnew=zeros(7,7);
    s=signal.data(1,:);     % since all transmissions have the same data
    Ntot=size(beams_ref,1);
    for k=1:Ntot            
        f=diag(beams_ref(k,:))*aBS;
        da=daUE.'*f;
        a=aUE.'*f;
        myGradient=zeros(signal.Ns,7);            
        g1=alphaLOS*s.*ddtauLOS;
        g2=alphaNLOS*s.*ddtauNLOS*a;
        g3=alphaNLOS*da*s.*dtauNLOS;
        g4=s.*dtauLOS;
        g5=1j*s.*dtauLOS;
        g6=a*s.*dtauNLOS;
        g7=1j*a*s.*dtauNLOS;
        myGradient=sqrt(signal.Es)*[g1.' g2.' g3.' g4.' g5.' g6.' g7.'];
        Jnew=Jnew+real(myGradient'*myGradient);
    end
    Jnew=Jnew/signal.sigma2;        
    Je=Jnew(1:3,1:3)-Jnew(1:3,4:end)*inv(Jnew(4:7,4:7))*Jnew(4:7,1:3); % EFIM of angles and delays
    J=Jnew;


    % Jacobian
    Gamma=zeros(3,3);
    eBS=(UE.location-BS.location)/norm(UE.location-BS.location);
    eRIS=(UE.location-RIS.location)/norm(UE.location-RIS.location);
    etheta=[0 -1; 1 0]*(UE.location-RIS.location)/(norm(UE.location-RIS.location)^2);
    c=signal.c;
    Gamma=[eBS/c eRIS/c etheta; 1/c 1/c 0];
    
    % FIM in position domain
    Jpos=Gamma*Je*Gamma';
    Jtmp=inv(Jpos);
    SPEB=Jtmp(1,1)+Jtmp(2,2);
    if (SPEB<0)
        SPEB=+inf;
    end
    PEB=sqrt(SPEB);    
    
