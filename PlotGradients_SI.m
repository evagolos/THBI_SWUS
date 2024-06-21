% Plot figs for the supplement


addpath('/Users/evagolos/Research/Misc/Colormaps/slanCM/');
addpath('/Users/evagolos/Research/ReceiverFunctions/THBI_from_Zach/SWUS/Visualization/github_repo/');


load data/Gradients_all.mat




%% Fig S5: Vs at NVGs, no interpretation

f1= figure;
ppos= get(f1,'Position');
set(f1,'Position',[ppos(1) ppos(2) ppos(3)*1.7 ppos(4)*3])

% Shallow NVG depth
subplot(3,2,1); hold on;
scatter(Lon(NVG1_Z>0),Lat(NVG1_Z>0),80,NVG1_Z(NVG1_Z>0),'filled');
scatter(Lon(gradual),Lat(gradual),80,'k','s','filled');
c= colorbar();
colormap(gca,flipud(slanCM(61)));
c.Label.String= 'Depth (km)';
caxis([60 250]);
plotRegions(gca);
xlim([-117.5 -104]); ylim([32 42.5]);
ylabel('Latitude (^o)');
set(gca,'FontSize',16);
title('a) Depth to shallow NVG');

% Shallow NVG Vs
subplot(3,2,2); hold on;
scatter(Lon(NVG1_Z_bottom>0),Lat(NVG1_Z_bottom>0),80,NVG1_Vs_bottom(NVG1_Z_bottom>0),'filled');
scatter(Lon(gradual),Lat(gradual),80,'s','k','filled');
c= colorbar();
cmap = getPyPlot_cMap('Spectral', 1000);
colormap(gca,cmap);
c.Label.String= 'Vs (km/s)';
plotRegions(gca);
xlim([-117.5 -104]); ylim([32 42.5]);
set(gca,'FontSize',16);
title('b) Vs below shallow NVG');

% Deep NVG depth
subplot(3,2,3); hold on;
isNZ2= NVG2_Z>0 & NVG2_Z<max(NVG1_Z_bottom)+10;
scatter(Lon(isNZ2),Lat(isNZ2),80,NVG2_Z(isNZ2),'filled');
c= colorbar();
colormap(gca,flipud(slanCM(61)));
c.Label.String= 'Depth (km)';
caxis([60 250]);
plotRegions(gca);
xlim([-117.5 -104]); ylim([32 42.5]);
xlabel('Longitude (^o)'); ylabel('Latitude (^o)');
set(gca,'FontSize',16);
title('c) Depth to deep NVG');

% Deep NVG Vs
subplot(3,2,4); hold on;
scatter(Lon(isNZ2),Lat(isNZ2),80,NVG2_Vs_bottom(isNZ2),'filled');
c= colorbar();
c.Label.String= 'Vs (km/s)';
cmap = getPyPlot_cMap('Spectral', 1000);
colormap(gca,cmap);
plotRegions(gca);
xlim([-117.5 -104]); ylim([32 42.5]);
xlabel('Longitude (^o)');
set(gca,'FontSize',16);
title('d) Vs below deep NVG');

% Reliability. Logical variable is1a denotes high Reliability.
subplot(3,2,5); hold on;
scatter(NVG1_Z_bottom(is1a),NVG1_Vs_bottom(is1a),50,[0.6350 0.0780 0.1840].*ones(sum(is1a),3),'filled');
scatter(NVG1_Z_bottom(~is1a),NVG1_Vs_bottom(~is1a),50,[0.9290 0.6940 0.1250].*ones(sum(~is1a),3),'filled');
title('e) Vs vs. Depth at base of first NVG');
xlabel('Depth to layer 1 base (km)'); ylabel('Vs (km/s)');
set(gca,'FontSize',16);
legend('High Reliability: Vs<=4.32 km/s','Low Reliability: Vs>4.32 km/s');

%% Figure S6: LAB depth + volcanic ages + MLD locations.

% Load volcanic information from NAVDAT.
AgesVolc= xlsread('data/navdat213821.xlsx','G5:G1271'); % "calculated" range
LatVolc= xlsread('data/navdat213821.xlsx','L5:L1271');
LonVolc= xlsread('data/navdat213821.xlsx','M5:M1271');
%AgeRange= AgesVolc(:,2)-AgesVolc(:,3);
%isAge= AgesVolc(:,1)~=0 & abs(AgeRange)<8;
isAge= AgesVolc>0 & AgesVolc<=80;
isYoung= AgesVolc~=0 & AgesVolc(:,1)<=1;

% Making the plot: There's gotta be a better way, but I just load this file to help index
% inversion points into a grid.
tbl= importdata('profile_pts_SWUS05.txt');
pts_lat = tbl(:,1)/10; pts_lon = tbl(:,2)/10;
len=length(pts_lat);
nlat= length(unique(pts_lat)); nlon= length(unique(pts_lon));
badlist= [533 534 561 562 563 564 565 566 589 590 591 592 593 594 595 596 597]; % No surface wave data.
LABgrid= zeros(nlon,nlat);

for ii=1:length(Lat)
    ind_list= find(pts_lon==Lon(ii) & pts_lat==Lat(ii));
    ilat= ceil(ind_list/nlon); ilon= mod(ind_list,nlon);
    if ilon==0; ilon=nlon; end
    LABgrid(ilon,ilat)= LAB_Z(ii);
end

figure; hold on;

% Put NVG1 depth into a regular grid to display as a smoothed surface.
ifill= LABgrid==0; LABgrid(ifill)= NaN;
LABgrid= fillmissing(LABgrid,"linear");
longrid= reshape(pts_lon,[nlon nlat]);
latgrid= reshape(pts_lat,[nlon nlat]);
% Now plot LAB surface.
apos= get(gca,'Position');
ax1= axes('Position',[apos(1)*0.7 0.2 apos(3)*0.9 apos(4)*0.89]); hold on;
[~,h]= contourf(longrid,latgrid,LABgrid,100);
set(h,'LineColor','none');
caxis(ax1,[60 250]);
colormap(gca,flipud(slanCM(61)));
plotRegions(ax1);
set(gca,'FontSize',16);
xlim([-117.5 -104]); ylim([32 42.5]);
% Volcano ages.
ax2= axes('Position',[apos(1)*0.7 0.2 apos(3)*0.9 apos(4)*0.89]); hold on;
scatter(ax2,LonVolc(isAge),LatVolc(isAge),30,AgesVolc(isAge),'d','filled','MarkerEdgeColor','k');
scatter(ax2,Lon(shallow_is_MLD),Lat(shallow_is_MLD),50,'k','^','filled');
ax2.XTick = [];  ax2.YTick = [];
linkaxes([ax1 ax2],'xy');
ax2.Visible = 'off';
colormap(ax2,'cool');
xlim(ax1,[-117.5 -104]); ylim([32 42.5]);
xlim(ax2,[-117.5 -104]); ylim([32 42.5]);
set(gca,'FontSize',14);
ylabel(ax1,'Latitude');
cb1= colorbar(ax1,'Location','southoutside','Position',[0.12 0.11 0.2 0.02]); 
cb1.Label.String= 'Depth (km)';
cb2= colorbar(ax2,'Location','southoutside','Position',[0.53 0.11 0.2 0.02]);
cb2.Label.String= 'Age (Ma)';
title(ax1,'Smooth LAB depth and Magmatic Age');


%% Fig S7: Moho depth and uncertainty
f4= figure;
ppos= get(f4,'Position');
set(f4,'Position',[ppos(1) ppos(2) ppos(3)*2 ppos(4)])

% Moho Depth
subplot(1,2,1); hold on;
scatter(Lon(Moho_Z>0),Lat(Moho_Z>0),80,Moho_Z(Moho_Z>0),'filled');
c= colorbar();
colormap(gca,flipud(slanCM(61)));
c.Label.String= 'km';
caxis([25 60]);
plotRegions(gca);
xlim([-117.5 -104]); ylim([32 42.5]);
set(gca,'FontSize',16);
xlabel('Longitude (^o)'); ylabel('Latitude (^o)');
title('c) Moho Depth');


% Moho uncertainty
subplot(1,2,2); hold on;
% Load models to get posterior distribution mu and sigma.
fresults= 'mdls_SWUS.mat';
load(fresults);
lenm= length(mdls.lat);
Moho_sigma= zeros(lenm,1); lon= zeros(lenm,1); lat= zeros(lenm,1);
% Load all profiles.
for ii=1:lenm
    % Find index for filling in grid.
    lon(ii)= mdls.lon(ii); lat(ii)= mdls.lat(ii);
    % Fill in each variable to plot.
    model= mdls.model{ii};
    Moho_sigma(ii)= model.Zd(2).std;
end
scatter(lon,lat,80,Moho_sigma,'filled');
c= colorbar();
colormap(gca,flipud(slanCM(61)));
c.Label.String= 'km';
plotRegions(gca);
xlim([-117.5 -104]); ylim([32 42.5]);
set(gca,'FontSize',16);
xlabel('Longitude (^o)'); ylabel('Latitude (^o)');
title('c) Moho Depth Standard Deviation');
