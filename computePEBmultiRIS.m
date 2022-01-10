function PEBtot=computePEBmultiRIS(BS,UE,signal,JLOS,JNLOS,locations,L)
% function PEBtot=computePEBmultiRIS(BS,UE,signal,JLOS,JNLOS,locations,L)
% (c) 2021, Henk Wymeersch, henkw@chalmers.se
% This is version 1.0 (Last edited: 2022-01-10)
% 
% Operation: combines FIM from several RIS to compute a PEB
% 
% inputs:
%   -BS: BS structure
%   -UE: UE structure
%   -signal: signal structure
%   -JLOS: LOS FIM (BS->UE)
%   -JNLOS: FIM for each RIS
%   -locations: 2D locations of each RIS
%   -L: number of RIS

% create the large overall FIM, including LOS path and all NLOS (RIS) paths
J=zeros(3+4*L,3+4*L);
J(1:3,1:3)=JLOS;                        %contains: [tauLOS, R{alphaLOS}, I{alphaLOS}]
ii=4;
for l=1:L
    J(ii:ii+3,ii:ii+3)=JNLOS(:,:,l);    %contains: [tauNLOS, theta, R{alphaNLOS}, I{alphaNLOS}]
    ii=ii+4;
end

% create the large Jacobian
eBS=(UE.location-BS.location)/norm(UE.location-BS.location);
c=signal.c;
Gamma=zeros(3+(L+1)*2,3+4*L);
Gamma(1:3,1)=[eBS/c; 1/c];
Gamma(4:5,2:3)=eye(2);
ii=4;
for l=1:L
    eRIS=(UE.location-locations(:,l))/norm(UE.location-locations(:,l));
    etheta=[0 -1; 1 0]*(UE.location-locations(:,l))/(norm(UE.location-locations(:,l))^2);
    Gamma(1:3,ii:ii+1)=[eRIS/c etheta; 1/c 0];
    Gamma((l-1)*2+6:(l-1)*2+7,ii+2:ii+3)=eye(2);
    ii=ii+4;
end
Jpos=Gamma*J*Gamma';
Jtmp=inv(Jpos);
SPEB=Jtmp(1,1)+Jtmp(2,2);
if (SPEB<0)
    SPEB=+inf;
end
PEBtot=sqrt(SPEB);   
