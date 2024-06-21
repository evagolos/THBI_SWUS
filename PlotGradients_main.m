% Plot main text Figures


addpath('/Users/evagolos/Research/Misc/Colormaps/slanCM/');
addpath('/Users/evagolos/Research/ReceiverFunctions/THBI_from_Zach/SWUS/Visualization/github_repo/');



load data/Gradients_all.mat





%% Figure 5: LAB and MLD maps
f0= figure;
ppos= get(f0,'Position');
set(f0,'Position',[ppos(1) ppos(2) ppos(3)*1.7 ppos(4)*1.8]);

% LAB depth
ax1= axes('Position',[0.08 0.55 0.4 0.4]); hold on;
scatter(Lon(is1a & ~gradual),Lat(is1a & ~gradual),80,NVG1_Z(is1a & ~gradual),'filled','MarkerEdgeColor','k','linewidth',1);
scatter(Lon(shallow_is_LAB),Lat(shallow_is_LAB),80,NVG1_Z(shallow_is_LAB),'filled');
scatter(Lon(deep_is_LAB & ~shallow_is_MLD),Lat(deep_is_LAB & ~shallow_is_MLD),80,NVG2_Z(deep_is_LAB & ~shallow_is_MLD),'filled');
scatter(Lon(shallow_is_MLD & deep_is_LAB),Lat(shallow_is_MLD & deep_is_LAB),80,NVG2_Z(shallow_is_MLD & deep_is_LAB),'^','filled');
scatter(Lon(gradual),Lat(gradual),80,'k','s','filled');
c= colorbar(); c.Label.String= 'Depth (km)';
caxis([60 250]);
colormap(gca,flipud(slanCM(61)));
plotRegions(gca);
xlim([-117.5 -104]); ylim([32 42.5]);
ylabel('Latitude (^o)');
set(gca,'FontSize',16);
title('a) Inferred Depth to LAB');

% MLD depth
ax2= axes('Position',[0.55 0.55 0.4 0.4]); hold on; hold on;
scatter(Lon(shallow_is_MLD),Lat(shallow_is_MLD),80,NVG1_Z(shallow_is_MLD),'^','filled');
c= colorbar(); c.Label.String= 'Depth (km)';
caxis([60 250]);
colormap(gca,flipud(slanCM(61)));
plotRegions(gca);
xlim([-117.5 -104]); ylim([32 42.5]);
%ylabel('Latitude (^o)'); xlabel('Longitude (^o)'); 
set(gca,'FontSize',16);
title('b) Inferred Depth to MLD');

% LAB Vs contrast
ax3= axes('Position',[0.08 0.08 0.4 0.4]); hold on; hold on;
scatter(Lon(is1a & ~gradual),Lat(is1a & ~gradual),80,LAB_amp(is1a & ~gradual),'filled','MarkerEdgeColor','k','linewidth',1);
scatter(Lon(shallow_is_LAB),Lat(shallow_is_LAB),80,LAB_amp(shallow_is_LAB),'filled');
scatter(Lon(deep_is_LAB & ~shallow_is_MLD),Lat(deep_is_LAB & ~shallow_is_MLD),80,LAB_amp(deep_is_LAB & ~shallow_is_MLD),'filled');
scatter(Lon(deep_is_LAB & shallow_is_MLD),Lat(deep_is_LAB & shallow_is_MLD),80,LAB_amp(deep_is_LAB & shallow_is_MLD),'^','filled');
scatter(Lon(gradual),Lat(gradual),80,'k','s','filled');
c= colorbar(); c.Label.String= '% dVs Contrast';
caxis([0 10]);
plotRegions(gca);
xlim([-117.5 -104]); ylim([32 42.5]);
xlabel('Longitude (^o)'); ylabel('Latitude (^o)');
set(gca,'FontSize',16);
title('c) LAB Vs Contrast');

% LAB width
ax4= axes('Position',[0.55 0.08 0.4 0.4]); hold on; hold on;
scatter(Lon(is1a & ~gradual),Lat(is1a & ~gradual),80,LAB_width(is1a & ~gradual),'filled','MarkerEdgeColor','k','linewidth',1);
scatter(Lon(shallow_is_LAB),Lat(shallow_is_LAB),80,LAB_width(shallow_is_LAB),'filled');
scatter(Lon(deep_is_LAB & ~shallow_is_MLD),Lat(deep_is_LAB & ~shallow_is_MLD),80,LAB_width(deep_is_LAB & ~shallow_is_MLD),'filled');
scatter(Lon(deep_is_LAB & shallow_is_MLD),Lat(deep_is_LAB & shallow_is_MLD),80,LAB_width(deep_is_LAB & shallow_is_MLD),'^','filled');
scatter(Lon(gradual),Lat(gradual),80,'k','s','filled');
c= colorbar(); c.Label.String= 'Depth (km)';
caxis([40 150]);
colormap(gca,flipud(slanCM(61)));
plotRegions(gca);
xlim([-117.5 -104]); ylim([32 42.5]);
xlabel('Longitude (^o)'); %ylabel('Latitude (^o)');
set(gca,'FontSize',16);
title('d) LAB Width (km)');



%% Figure 6: LAB width vs. depth trend

figure; hold on;
a1= scatter(LAB_Z(is1a & ~gradual),LAB_width(is1a & ~gradual),40,NVG1_Vs(is1a & ~gradual),'filled','MarkerEdgeColor','k','linewidth',1.2);
a2= scatter(LAB_Z(shallow_is_LAB),LAB_width(shallow_is_LAB),50,NVG1_Vs(shallow_is_LAB),'filled','MarkerEdgeColor',[0.5 0.5 0.5],'linewidth',0.5);
a3= scatter(LAB_Z(deep_is_LAB & ~shallow_is_MLD),LAB_width(deep_is_LAB & ~shallow_is_MLD),50,NVG2_Vs(deep_is_LAB & ~shallow_is_MLD),'filled','MarkerEdgeColor',[0.5 0.5 0.5],'linewidth',0.5);
a4= scatter(LAB_Z(deep_is_LAB & shallow_is_MLD),LAB_width(deep_is_LAB & shallow_is_MLD),50,NVG2_Vs(deep_is_LAB & shallow_is_MLD),'filled','^','MarkerEdgeColor',[0.5 0.5 0.5],'linewidth',0.5);
c= colorbar();
cmap = getPyPlot_cMap('Spectral', 1000);
colormap(gca,cmap);
caxis([4.18 4.7]);
c.Label.String= 'Vs (km/s)';
legend([a1 a4],'Only LAB','LAB with MLD');
xlabel('LAB Depth (km)'); ylabel('LAB gradient width (km)');
title('LAB Width vs. Depth');
set(gca,'FontSize',16);





%% Figure 7: PVG points below known LAB.

f3= figure;
ppos= get(f3,'Position');
set(f3,'Position',[ppos(1) ppos(2) ppos(3)*3 ppos(4)*2]);

% Histogram of PVG depths.
subplot(2,3,1); hold on;
histogram(PVG_Z(isPVG_amp),15);
xlabel('PVG Depth (km)'); ylabel('Number');
set(gca,'FontSize',16);
title('a) PVG Depth Distribution');

% Non-convincing plot of contrast vs. depth
subplot(2,3,2); hold on;
scatter(PVG_Z(isPVG_amp & PVG_amp>0),100*PVG_amp(isPVG_amp & PVG_amp>0),'filled');
xlabel('Depth to PVG center'); ylabel('PVG Vs Contrast');
set(gca,'FontSize',16);
title('b) Vs Contrast and Depth of PVG');

% Slightly convincing plot of contrast vs. Vs above PVG
subplot(2,3,3); hold on;
scatter(PVG_Vs_top(isPVG_amp & PVG_amp>0),100*PVG_amp(isPVG_amp & PVG_amp>0),'filled');
% Fit a line showing trend.
mdl= fitlm(PVG_Vs_top(isPVG_amp & PVG_amp>0),100*PVG_amp(isPVG_amp & PVG_amp>0));
xest= 3.9:0.1:4.6;
PVG_amp_est= mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2)*xest;
plot(xest,PVG_amp_est,'linewidth',1.2,'Color','k');
xlabel('Vs above PVG'); ylabel('PVG Vs Contrast');
set(gca,'FontSize',16);
title('c) Vs Contrast and Vs Above PVG');

% Depth to center of PVG
subplot(2,3,4); hold on;
scatter(Lon(isPVG_amp),Lat(isPVG_amp),80,PVG_Z(isPVG_amp),'filled');
scatter(Lon(isPVG_amp & is1a),Lat(isPVG_amp & is1a),120,PVG_Z(isPVG_amp & is1a),'filled','MarkerEdgeColor','k','linewidth',1.2);
c= colorbar();
c.Label.String= 'Depth (km)';
colormap(gca,flipud(slanCM(61)));
plotRegions(gca);
xlim([-117.5 -104]); ylim([32 42.5]);
set(gca,'FontSize',16);
xlabel('Longitude (^o)'); ylabel('Latitude (^o)');
title('d) Depth to PVG Center');

% Vs at the top of PVG
subplot(2,3,5); hold on;
scatter(Lon(isPVG_amp),Lat(isPVG_amp),80,PVG_Vs_bottom(isPVG_amp),'filled');
scatter(Lon(isPVG_amp & is1a),Lat(isPVG_amp & is1a),120,PVG_Vs_top(isPVG_amp & is1a),'filled','MarkerEdgeColor','k','linewidth',1.2);
c= colorbar();
c.Label.String= 'Vs (km/s)';
caxis([4 4.35]);
cmap = getPyPlot_cMap('Spectral',1000);
colormap(gca,cmap);
plotRegions(gca);
xlim([-117.5 -104]); ylim([32 42.5]);
set(gca,'FontSize',16);
xlabel('Longitude (^o)'); ylabel('Latitude (^o)');
title('e) Vs at top of PVG');

% Vs contrast over PVG.
subplot(2,3,6); hold on;
scatter(Lon(isPVG_amp),Lat(isPVG_amp),80,PVG_amp(isPVG_amp),'filled');
scatter(Lon(isPVG_amp & is1a),Lat(isPVG_amp & is1a),120,100*PVG_amp(isPVG_amp & is1a),'filled','MarkerEdgeColor','k','linewidth',1.2);
c= colorbar();
c.Label.String= '% Vs Contrast';
plotRegions(gca);
xlim([-117.5 -104]); ylim([32 42.5]);
set(gca,'FontSize',16);
xlabel('Longitude (^o)'); ylabel('Latitude (^o)');
title('f) PVG Vs Contrast');