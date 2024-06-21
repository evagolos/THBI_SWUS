% Plot data error from inversion hyperparameters, and model uncertainty for
% several parameters.

% Modified from plotting scripts of Brunsvik et al. (2024), GRL.

close all
clear


addpath('../../Visualization/github_repo/');
addpath('/Users/evagolos/Research/Misc/Colormaps/slanCM/');
addpath('/Users/evagolos/Research/ReceiverFunctions/lib/m_map/');
addpath('functions');

% Packages needed for plotting:
% slanCM: https://www.mathworks.com/matlabcentral/fileexchange/120088-200-colormap
% PyColormap4Matlab: https://www.mathworks.com/matlabcentral/fileexchange/68239-pycolormap4matlab
% m_map: https://www.eoas.ubc.ca/~rich/map.html

fresults= 'mdls_SWUS.mat';
load(fresults);
lenm= length(mdls.lat);

% For making grid.
tbl= importdata('profile_pts_SWUS05.txt');
pts_lat = tbl(:,1)/10; pts_lon = tbl(:,2)/10;
len=length(pts_lat);
nlat= length(unique(pts_lat)); nlon= length(unique(pts_lon));
badlist= [533 534 561 562 563 564 565 566 589 590 591 592 593 594 595 596 597]; 

Z_toUse= [60 100 180];

% Locations for profiles in Figure 2.
locations= [-111.5 42; -109.5 39.5; -114.5 39.0; -111.0 37.0; -109.5 33.0; -105.0 34.5];

Vs= zeros(nlon,nlat,3);
Vs_sigma60= zeros(lenm,1); Vs_sigma100= zeros(lenm,1); Vs_sigma180= zeros(lenm,1);
lon= zeros(lenm,1); lat= zeros(lenm,1);

% Load all profiles.
for ii=1:lenm
    
    % Find index for filling in grid.
    lon(ii)= mdls.lon(ii); lat(ii)= mdls.lat(ii);
    ind_list= find(pts_lon==lon(ii) & pts_lat==lat(ii));
    ilat= ceil(ind_list/nlon); ilon= mod(ind_list,nlon);
    if ilon==0; ilon=nlon; end
    
    % Fill in each variable to plot.
    model= mdls.model{ii};
    Vs_sigma_prof= (model.Vssig2(:,2) - model.Vssig2(:,1))/4;
    Vs_sigma(ii,:)= interp1(model.Z,Vs_sigma_prof,Z_toUse,'linear');
    %Moho_sigma(ii)= model.Zd(2).std;
    Vs_sigma_prof= (model.Vssig2(:,2) - model.Vssig2(:,1))/4;
    Vs_sigma(ii,:)= interp1(model.Z,Vs_sigma_prof,Z_toUse,'linear');
    Vs_sigma60(ii)= Vs_sigma(ii,1);
    Vs_sigma100(ii)= Vs_sigma(ii,2);
    Vs_sigma180(ii)= Vs_sigma(ii,3);
    
    % Data errors. Fill these too while we're at it!
    sigma_CCP(ii)= model.hyperparms.sig_RF_Sp_CCP.mu_log10;
    sigma_SW(ii)= model.hyperparms.sig_SW_Ray_phV.mu_log10;
    
    % Vs
    %Vs_prof= (model.Vssig2(:,2) - model.Vssig2(:,1))/4;
    %
    if ~ismember(ind_list,badlist)
    %Vs_prof= nodel.; 
        Vs(ilon,ilat,1)= interp1(model.Z,model.Vsav,Z_toUse(1),'linear');
        Vs(ilon,ilat,2)= interp1(model.Z,model.Vsav,Z_toUse(2),'linear');
        Vs(ilon,ilat,3)= interp1(model.Z,model.Vsav,Z_toUse(3),'linear');
    else
        Vs(ilon,ilat,:)= NaN;
    end
    
end

% % For profiles along a given latitude line.
profs= [41 39 37 35 33];
nprofs= length(profs);
proflabel= {'A','B','C','D','E'};

% Establish grid -- need to populate this to make maps (panels a-c)
longrid= reshape(pts_lon,[nlon nlat]);
latgrid= reshape(pts_lat,[nlon nlat]);

%% Plot
figure; clf; set(gcf,'pos', [87 856 692*1.7 476*1.7]); 

n_plots= 9;
for ifig = 1:n_plots

if ifig<=3    
    var= Vs(:,:,ifig); cmap = getPyPlot_cMap('Spectral', 1000);
    cbarstr= 'Vs (km/s)';
    if ifig==1
        titlestr=sprintf('a) Vs at %.0f km',Z_toUse(ifig));
        axes('Position',[0.05 0.57 0.275 0.36]);
        ylabel('Latitude (^o)');
        caxis([4.2 4.8]);
    elseif ifig==2
        titlestr=sprintf('b) Vs at %.0f km',Z_toUse(ifig));
        axes('Position',[0.37 0.57 0.275 0.36])
        caxis([4.2 4.8]);
    elseif ifig==3
        titlestr=sprintf('c) Vs at %.0f km',Z_toUse(ifig));
        axes('Position',[0.7 0.57 0.275 0.36])
        caxis([4.2 4.8]);
    end
elseif ifig==4
    var= Vs_sigma60; cmap = getPyPlot_cMap('Spectral', 1000);
    titlestr= 'd) Vs stdev at 60 km';
    cbarstr= 'Vs (km/s)';
    ax= axes('Position',[0.05 0.1 0.275 0.36]);
    xlabel('                                Longitude (^o)'); 
    ylabel('Latitude (^o)');
elseif ifig==5
    var= Vs_sigma100; cmap = getPyPlot_cMap('Spectral', 1000);
    titlestr= 'e) Vs stdev at 100 km';
    cbarstr= 'Vs (km/s)';
    ax= axes('Position',[0.37 0.1 0.275 0.36]);
    xlabel('                                Longitude (^o)');
else
    var= Vs_sigma180; cmap = getPyPlot_cMap('Spectral', 1000);
    titlestr= 'f) Vs stdev at 180 km';
    cbarstr= 'Vs (km/s)';
    ax= axes('Position',[0.7 0.1 0.275 0.36]);
    xlabel('                                Longitude (^o)');
end

hold on; box on; set(gca, 'LineWidth', 1.5);


if ifig<=3
%    Vs_smooth= imgaussfilt(var,0.2);
    if ifig==1
        Vs_smooth= smoothdata(var,'SmoothingFactor',0.01);   
    else
        Vs_smooth= smoothdata(var,'SmoothingFactor',0.001);
    end
    
    [~,h]= contourf(gca,longrid,latgrid,Vs_smooth,3.9:0.005:4.8); set(h,'LineColor','none');
else
    scatter(lon,lat,80,var,'filled');
end

colormap(gca,cmap); 

c= colorbar('Location','southoutside'); 
c.Position(2)= c.Position(2)-0.075;
c.Position(3)= c.Position(3)*0.5;
c.Label.String= cbarstr;
c.Label.FontSize= 12;
set(gca,'FontSize',12);
title(titlestr,'FontSize',14);

plotRegions(gca);

if ifig==1        
        % Add profile lines and markers
        for iline=1:length(profs)
        plot([-117.5 -104],[profs(iline) profs(iline)],'-o','color','k','linewidth',2,'markerfacecolor','k','markersize',8);
        text(-117.5,profs(iline)+0.4,[proflabel{iline}],'FontSize',14);
        text(-104.5,profs(iline)+0.4,[proflabel{iline},''''],'FontSize',14);
        end
        
        for iloc= 1:length(locations)
            scatter(locations(iloc,1),locations(iloc,2),160,'g','filled','Marker','d'); 
            text(locations(iloc,1)+0.15,locations(iloc,2)+0.3,num2str(iloc),'FontSize',14);
        end
end
    
xlim([-117.5 -104]); ylim([32 42.5]); 

if ifig>3 
% Plot PDF of stdevs at all locations, overlapping the colorbar.  
ax_temp = axes('Position',c.Position); hold on; cla; 
[pdf_y,pdf_x] = ksdensity(var); 
plot(pdf_x, pdf_y, 'k', 'LineWidth',3); 
set(ax_temp, 'Color' , 'None', ...
    'XLim', c.Limits, ...
    'XTickLabel', [], 'YTickLabel', [] ); 
ylabel('P(m)'); 
disp('One plot done')
end
 
 
end
 
 

%% Figure S4: Data error for CCP stack and surface waves.
figure(2); clf; set(gcf,'pos', [87 856 692*1.2 476*0.88]); 
n_plots= 2;

for ifig=1:n_plots
if ifig==1
    var= sigma_CCP; cmap= 'parula';
    titlestr= 'a) CCP stack error';
    cbarstr= '\sigma';
    ax= axes('Position',[0.05 0.18 0.4 0.75]);
    ylabel('Latitude (^o)'); 
    xlabel('                                Longitude (^o)'); % So we can read the label.
elseif ifig==2
    var= sigma_SW; cmap= 'parula';
    titlestr= 'b) Surface wave error';
    cbarstr= '\sigma';
    ax= axes('Position',[0.55 0.18 0.4 0.75]);
    xlabel('                                Longitude (^o)');
    
end
 

hold on; box on; set(gca, 'LineWidth', 1.5);


scatter(lon,lat,80,var,'filled');


if ifig==3
    colormap(gca,flipud(slanCM(61)));
else
    colormap(gca,cmap); 
end
c= colorbar('Location','southoutside'); 
c.Position(2)= c.Position(2)-0.15;
c.Position(3)= c.Position(3)*0.5;
c.Label.String= cbarstr;
c.Label.FontSize= 12;
set(gca,'FontSize',12);
title(titlestr,'FontSize',14);

% Add tectonic boundaries.
plotRegions(gca);
    
xlim([-117.5 -104]); ylim([32 42.5]); 
 
% Plot PDF of sigma values at all locations, overlapping the colorbar. 
ax_temp = axes('Position',c.Position); hold on; cla; 
[pdf_y,pdf_x] = ksdensity(var); 
plot(pdf_x, pdf_y, 'k', 'LineWidth',3); 
set(ax_temp, 'Color' , 'None', ...
    'XLim', c.Limits, ...
    'XTickLabel', [], 'YTickLabel', [] ); 
ylabel('P(m)'); 
disp('One plot done');



end