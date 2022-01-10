function []=plotPEBContour(signal,RIS,BS,RISstate,x_grid,y_grid)
% function []=plotPEBContour(signal,RIS,BS,RISstate,x_grid,y_grid)
% (c) 2021, Henk Wymeersch, henkw@chalmers.se
% This is version 1.0 (Last edited: 2022-01-10)
% 
% Operation: Generates a PEB contour plot for fixed BS location and RIS located in
% RISstate over a grid x_grid x y_grid
% 
% inputs: 
%   -signal: signal structure
%   -RIS: RIS structure
%   -BS: BS structure
%   -RISstate: RISstate is a 3 x N_RIS matrix, with first row the X-coordinates of the RIS, second row the Y-coordinates of the RIS, and third row the orientation of the RIS
%   -x_grid: X-coordinates of the PEB evaluation
%   -y_grid: Y-coordinates of the PEB evaluation

L=size(RISstate,2);         % number of RIS
PEBcomb=zeros(length(x_grid),length(y_grid));
hw = waitbar(0,'Running Offline PEB Optimization ...');
for ix=1:length(x_grid)      
    waitbar(ix/length(x_grid),hw);
    for iy=1:length(y_grid)        
        UE.location=[x_grid(ix); y_grid(iy)];        
        % compute the PEB
        JLOS=zeros(3,3);
        JNLOS=zeros(4,4,L);
        for is=1:L    
            RIS.location=RISstate(1:2,is); 
            RIS.kappa=RISstate(3,is);
            locations(:,is)=RIS.location;
            orientations(:,is)=RIS.kappa;
            beams_ref=exp(1j*rand(signal.T,RIS.M)*2*pi);        % generate T random beams   
            [~,J]=computePEBSPM(BS,UE,RIS,signal,beams_ref);    % order [tauLOS, tauNLOS, theta, R{alphaLOS}, I{alphaLOS}, R{alphaNLOS}, I{alphaNLOS}]
            JLOS=J([1 4 5],[1 4 5]);
            JNLOS(:,:,is)=J([2 3 6 7],[2,3,6,7]);           
        end            
        % now combine information from all the RIS to compute the PEB
        PEBcomb(ix,iy)=computePEBmultiRIS(BS,UE,signal,JLOS,JNLOS,locations,L);        
    end   
end
close(hw)
tmp=PEBcomb(:);
drawArrow = @(x,y) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0.5,'Color','w','LineWidth',2,'AutoScale','off');
figure
e=[1;0];
f3=contourf(x_grid,y_grid,10*log10(min(PEBcomb,20))','edgecolor','none');
hold on
plot(BS.location(1),BS.location(2),'rs')
tt=text(BS.location(1)-0.5,BS.location(2)+0.3,'BS');
set(tt,'Interpreter','latex','FontSize',12);
for is=1:L    
    x1=locations(1,is);
    y1=locations(2,is);
    R=[cos(orientations(is)) -sin(orientations(is)); sin(orientations(is)) cos(orientations(is))];
    u=R'*2*e;
    x2=x1+u(1);
    y2=y1+u(2);
    plot(locations(1,is),locations(2,is),'b*')
    tt=text(locations(1,is)-u(1)/4-0.3,locations(2,is)-u(2)/4,'RIS');
    set(tt,'Interpreter','latex','FontSize',12);
    drawArrow([x1 x2],[y1 y2])    
end
rectangle('Position',[-5 0 10 10],'Edgecolor','b','LineWidth',2)
colorbar
xl=xlabel('$X$ coordinate [m]');
yl=ylabel('$Y$ coordinate [m]');
tl=title('PEB [dB-meter]');
set(xl,'Interpreter','latex','FontSize',12);
set(yl,'Interpreter','latex','FontSize',12);
set(tl,'Interpreter','latex','FontSize',12);
set(gcf, 'Color', 'w');
set(gca,'fontsize',16);


