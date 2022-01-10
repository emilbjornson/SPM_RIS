function onlinedesign(BS,UE,RIS,signal)
% function onlinedesign(BS,UE,RIS,signal)
% (c) 2021, Henk Wymeersch, henkw@chalmers.se
% This is version 1.0 (Last edited: 2022-01-10)
% 
% Evaluates the PEB for different RIS phase profiles
% inputs: 
%   -BS: BS structure
%   -UE: UE structure
%   -RIS: RIS structure
%   -signal: signal structure

    % generate response vectors and derivatives
    [~,aUE,daUE,~]=computeRISChannelSPM(UE.location,RIS,signal);                  % channel from RIS to UE
    [~,aBS,~,~]=computeRISChannelSPM(BS.location,RIS,signal);                     % channel from BS to RIS
   
    % now generate the combined response vector and the derivatives, which will be
    % used as RIS phase profiles
    b=aUE.*aBS;             % combined response BS->RIS and RIS->UE
    db=daUE.*aBS;           % derivative with respect to AOD
    db=db*norm(b)/norm(db); % normalize the derivative beam

    % under the optimal design we only use the following 2 beam
    beams_refOD(1,:)=conj(b);   % directional vector
    beams_refOD(2,:)=conj(db);  % derivative vector

    Ntot=256;   % let's use Ntot transmission
    % now sweep over the time sharing factor and determine the PEB
    hw = waitbar(0,'Running Online PEB Optimization ...');
    for n1=0:Ntot           % number of times we use directional beam            
        waitbar(n1/Ntot,hw);
        % compute the FIM and PEB                
        beamsA=repmat(beams_refOD(1,:),n1,1);
        beamsB=repmat(beams_refOD(2,:),Ntot-n1,1);        
        beams=[beamsA; beamsB];
        [PEB(n1+1),Jnew]=computePEBSPM(BS,UE,RIS,signal,beams); 
        Je=Jnew(1:3,1:3)-Jnew(1:3,4:end)*inv(Jnew(4:7,4:7))*Jnew(4:7,1:3); % EFIM of angles and delays 
        Se=inv(Je);
        CRB(n1+1,1)=sqrt(Se(1,1));
        CRB(n1+1,2)=sqrt(Se(2,2))*signal.c;
        CRB(n1+1,3)=sqrt(Se(3,3));  
    end
    close(hw)
    % compare with PEB of random beams
    beams_random=exp(1j*rand(Ntot,RIS.M)*2*pi);        % generate random beams   
    [PEB_random,Jnew]=computePEBSPM(BS,UE,RIS,signal,beams_random); 
    Je=Jnew(1:3,1:3)-Jnew(1:3,4:end)*inv(Jnew(4:7,4:7))*Jnew(4:7,1:3); % EFIM of angles and delays           
    Se=inv(Je);    
    CRB_random(1)=sqrt(Se(1,1));
    CRB_random(2)=sqrt(Se(2,2))*signal.c;
    CRB_random(3)=sqrt(Se(3,3));
   

    % now visualize the results
    kk=(0:Ntot)/Ntot;    
    figure;
    semilogy(kk,PEB,'k-',kk,ones(1,Ntot+1)*PEB_random,'k--','Linewidth',2)
    hold on
    semilogy(kk,CRB(1:end,2),'b-',kk,ones(1,Ntot+1)*CRB_random(2),'b--','Linewidth',2);
    semilogy(kk,CRB(1:end,3),'r-',kk,ones(1,Ntot+1)*CRB_random(3),'r--','Linewidth',2);
    xl=xlabel('$T_1/T$');
    yl=ylabel('accuracy');
    ll=legend('PEB [m] designed','PEB [m] random','RIS-TOA error [m] designed','RIS-TOA error [m] random','RIS-AOD error [rad] designed','RIS-AOD error [rad] random');
    set(xl,'Interpreter','latex','FontSize',12);
    set(yl,'Interpreter','latex','FontSize',12);
    set(ll,'Interpreter','latex','FontSize',16,'Location','SouthEast');
    grid
    axis([0 1 1e-5 1e-1])
    set(gca,'fontsize',16);
    set(gcf, 'Color', 'w');

 
   
    