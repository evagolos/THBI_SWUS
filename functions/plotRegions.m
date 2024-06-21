function [] = plotRegions(ax)
% Plots geographic boundaries, EG 2021. 
% ax = enter axis label to plot on (if none specified, use gca)

%hold on;    
BNR= importdata('data/GeoBounds/BNR.txt');
Col= importdata('data/GeoBounds/Colorado.txt');
RGR= importdata('data/GeoBounds/RioGrande.txt');
%SNB= importdata('data/GeoBounds/Sierra.txt');
WY= importdata('data/GeoBounds/Wyoming.txt');

S=shaperead('usastatelo','UseGeoCoords',true);
for i=3:51; plot(S(i).Lon,S(i).Lat,'k'); end
plot(ax,Col(:,1),Col(:,2),'k','linewidth',2);
    text(ax,-112,37.5,'CP','FontWeight','bold','FontSize',14); % Can move text around.
plot(ax,RGR(:,1),RGR(:,2),'k','linewidth',2);
    text(ax,-105.8,32.5,'RG','Fontweight','bold','FontSize',12);
    %plot([-106 -105.5],[31 32],'k','linewidth',2);
plot(ax,BNR(:,1),BNR(:,2),'k','linewidth',2);
    text(ax,-116,40.1,'BNR','FontWeight','bold','FontSize',14);
% plot(ax,SNB(:,1),SNB(:,2),'k','linewidth',2);
%     text(ax,-120.8,38.2,'SN','FontWeight','bold','FontSize',12);
plot(ax,WY(:,1),WY(:,2),'k','linewidth',2);
    text(ax,-110,42,'WY','FontWeight','bold','FontSize',14);

 xlim([-130 -95]); ylim([30 50]); % Can change limits of plot

end
