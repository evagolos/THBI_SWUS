% Get gradient information from inverted models.

addpath('functions/');

fresults= 'mdls_SWUS.mat';
load(fresults);


% A whole bunch of variables!
len= length(mdls.model);
Vs_Z1= zeros(len,1); Vs_Z2= zeros(len,1); Vs_Z3= zeros(len,1);
isGradual= zeros(len,1);

Moho_Z= zeros(len,1); Moho_Vs= zeros(len,1); 
NVG1_Z= zeros(len,1); NVG1_Vs= zeros(len,1); % NVG1 refers to shallower NVG, or only NVG if one.
NVG1_Z_top= zeros(len,1); NVG1_Vs_top= zeros(len,1);
NVG1_Z_bottom= zeros(len,1); NVG1_Vs_bottom= zeros(len,1); 
NVG2_Z= zeros(len,1); NVG2_Vs= zeros(len,1); % NVG2 is the deeper NVG. Info is zero if it doesn't exist.
NVG2_Z_top= zeros(len,1); NVG2_Vs_top= zeros(len,1); 
NVG2_Z_bottom= zeros(len,1); NVG2_Vs_bottom= zeros(len,1);
NVG3_Z= zeros(len,1); NVG3_Vs= zeros(len,1); % A few points have three!

PVG_Vs= zeros(len,1);
PVG_Z= zeros(len,1); PVG_amp= zeros(len,1); % If the last positive peak exists, so does PVG.
PVG_Z_bottom= zeros(len,1); % Record info about last positive peak, if exists.
PVG_Vs_bottom= zeros(len,1);
PVG_Z_top= zeros(len,1); PVG_Vs_top= zeros(len,1);


category= cell(len,1);

% Loop through all models.
for ii=1:len
    
    % Load that model. Have to deal with each at a time; Z is different for each model because of
    % transdimensional nature of inversion.
    final_model= mdls.model{ii};
    Z= final_model.Z; Vs= final_model.Vsav;
    
    % Now we go through and select depth, Vs, and width of gradients
    % The Moho is the easiest case -- this is its own model parameter. Just
    % need to find Vs at that depth.
    zmoh= final_model.Zd(2).mu;
    Moho_Z(ii)= zmoh; imoh= find(Z==final_model.Zd(2).mu);
    if isempty(imoh)
        Zdiff= abs(Z-final_model.Zd(2).mu); imoh= find(Zdiff==min(Zdiff)); imoh= imoh(1);
    end
    Moho_Vs(ii)= Vs(imoh+2);
    icrust= find(Vs<Moho_Vs(ii)*0.95 & Z<Moho_Z(ii));
    
    % Now for NVGs.
    % First, find peaks in Z profile, both positive (pklocs) and negative (Npklocs)
    [~,pklocs,pw,prom]= findpeaks(Vs,'MinPeakProminence',0.001);
    Z_p= Z(pklocs); Vs_p= Vs(pklocs);
    [~,Npklocs,npw,nprom]= findpeaks(-Vs,'MinPeakProminence',0.001);
    Z_np= Z(Npklocs); Vs_np= Vs(Npklocs);
    
    % Remove low-prominence/erratic/too-deep positive peaks.
    removeMeP= pw.*prom<0.0155 | (pw.*prom<0.5 & Z_p>250);
    if ii==108; removeMeP(2)= true; end
    Z_p(removeMeP)= []; Vs_p(removeMeP)= []; pklocs(removeMeP)= [];
    
    % Remove low-prominence/erratic negative peaks. Stricter than positive
    % peak criteria.
    removeMe= npw.*nprom<0.61 & nprom <0.025;
    
    % Special case: don't ignore this peak.
    if ii==188; removeMe(end)= false; end
    Z_np(removeMe)= []; Vs_np(removeMe)= []; Npklocs(removeMe)= [];
    
% Now, determine whether there's one or two NVGs.
% No negative velocity peaks lower than 10 km below Moho.
if isempty(Z_np(Z_np>zmoh+10))
    % Gradual LAB case.
    [nvg_z,nvg_w,nvg_a] = model_NVG_info(final_model); % This function provides an estimate
    
    NVG1_Z_bottom(ii)= nvg_z+nvg_w/2;
    iz= find(abs(Z-NVG1_Z_bottom(ii))==min(abs(Z-NVG1_Z_bottom(ii)))); iz= iz(1);
    NVG1_Vs_bottom(ii)= Vs(iz);
    
    NVG1_Z_top(ii)= nvg_z-nvg_w/2;
    iz= find(abs(Z-NVG1_Z_top(ii))==min(abs(Z-NVG1_Z_top(ii)))); iz= iz(1);
    NVG1_Vs_top(ii)= Vs(iz);
    
    % Center of gradient
    NVG1_Z(ii)= (NVG1_Z_top(ii) + NVG1_Z_bottom(ii))/2; 
    iz= find(abs(Z-NVG1_Z(ii))==min(abs(Z-NVG1_Z(ii)))); iz= iz(1);
    NVG1_Vs(ii)= Vs(iz);
    
    isGradual(ii)= 1;
    
    
% Positive peak below last negative peak.
elseif Z_p(end) > Z_np(end)

    gVs = gradient(Vs,Z);
    

    if length(Z_np(Z_np>zmoh+10)) == 1 % Two NVGs: one in peaks, one below last positive peak.
        % "Base of NVG" info. Depth below last NVG is halfway between last positive peak and
        % model bottom.
        NVG1_Z_bottom(ii)= Z_np(end);
        NVG1_Vs_bottom(ii)= Vs(Npklocs(end));
        
        % Find next positive peak above that NVG.
        ip= find(Z_p < Z_np(end)); ip= ip(end);
        NVG1_Z_top(ii)= Z_p(ip);
        NVG1_Vs_top(ii)= Vs_p(ip);
        
        % Now info on center of gradients.
        NVG1_Z(ii)= (Z_np(end)+Z_p(ip))/2;
        iz= find(abs(Z-NVG1_Z(ii))==min(abs(Z-NVG1_Z(ii)))); iz= iz(1);
        NVG1_Vs(ii)= Vs(iz);
        
        % Second gradient.
        NVG2_Z_bottom(ii)= Z(end); NVG2_Vs_bottom(ii)= Vs(end); 
        NVG2_Z_top(ii)= Z_p(end); NVG2_Vs_top(ii)= Vs_p(end);

        % Center of second NVG: use same info as "bottom"
        NVG2_Z(ii)= (Z_p(end)+300)/2;
        iz= find(abs(Z-NVG2_Z(ii))==min(abs(Z-NVG2_Z(ii)))); iz= iz(1);
        NVG2_Vs(ii)= Vs(iz);
        
        % Bottom PVG.
        ip= find(Z_p > Z_np(end)); ip= ip(1);
        PVG_Z(ii)= (Z_p(ip)+NVG1_Z_bottom(ii))/2;
        PVG_Z_top(ii)= NVG1_Z_bottom(ii); PVG_Vs_top(ii)= NVG1_Vs_bottom(ii);
        PVG_Z_bottom(ii)= Z_p(ip); PVG_Vs_bottom(ii)= Vs_p(ip);
        
        % A few special cases.
        % One NVG but a prominent PVG below.
        if ii==254 
            NVG2_Z_bottom(ii)= 0; NVG2_Vs_bottom(ii)= 0;
            NVG2_Vs(ii)= 0; NVG2_Z(ii)= 0;
            PVG_Z(ii)= (220+182)/2; 
            PVG_Z_top(ii)= 182; PVG_Z_bottom(ii)= 220;
            iz= find(abs(Z-182)==min(abs(Z-182)));
            PVG_Vs_top(ii)= Vs(iz); 
            iz= find(abs(Z-220)==min(abs(Z-220)));
            PVG_Vs_bottom(ii)= Vs(iz);
        elseif ii==81
            NVG2_Z_bottom(ii)= 0; NVG2_Vs_bottom(ii)= 0;
            NVG2_Vs(ii)= 0; NVG2_Z(ii)= 0;
            PVG_Z_top(ii)= 274; % It's really quite confined to lower model here.
            PVG_Z_bottom(ii)= 300;
            iz= find(abs(Z-PVG_Z_top(ii))==min(abs(Z-PVG_Z_top(ii))));
            PVG_Vs_top(ii)= Vs(iz); 
            iz= find(abs(Z-PVG_Z_bottom(ii))==min(abs(Z-PVG_Z_bottom(ii))));
            PVG_Vs_bottom(ii)= Vs(iz);
        end
            
        iz= find(abs(Z-PVG_Z(ii))==min(abs(Z-PVG_Z(ii)))); iz= iz(1);
        PVG_Vs(ii)= Vs(iz);
        PVG_amp(ii)= (PVG_Vs(ii) - NVG1_Vs_bottom(ii))/mean([PVG_Vs(ii) NVG1_Vs_bottom(ii)]);
        
    else % It could happen!
        % NVG1_Z_bottom happens at last Z_np. NVG2_Z_bottom happens below last positive peak.
        NVG1_Z_top(ii)= Z_p(1); NVG1_Vs_top(ii)= Vs_p(1);
        NVG1_Z_bottom(ii)= Z_np(end-1); NVG1_Vs_bottom(ii)= Vs(Npklocs(end-1));
        NVG2_Z_top(ii)= Z_p(2); NVG2_Vs_top(ii)= Vs_p(2);
        NVG2_Z_bottom(ii)= Z_np(end); NVG2_Vs_bottom(ii)= Vs(Npklocs(end));
        NZ3(ii)= (Z_p(end)+300)/2;
        iz= find(abs(Z-NZ3(ii))==min(abs(Z-NZ3(ii)))); iz= iz(1);
        Vs_Z3(ii)= Vs(iz);
        
        % Center of gradients.
        NVG1_Z(ii)= (Z_np(end-1)+Z_p(1))/2;
        iz= find(abs(Z-NVG1_Z(ii))==min(abs(Z-NVG1_Z(ii)))); iz= iz(1);
        NVG1_Vs(ii)= Vs(iz);
        NVG2_Z(ii)= (Z_np(end)+Z_p(2))/2;
        iz= find(abs(Z-NVG2_Z(ii))==min(abs(Z-NVG2_Z(ii)))); iz= iz(1);
        NVG2_Vs(ii)= Vs(iz);
        NVG3_Z(ii)= NZ3(ii); NVG3_Vs= Vs_Z3(ii);
        
        PVG_Z_top(ii)= NVG2_Z_bottom(ii); PVG_Vs_top(ii)= NVG2_Vs_bottom(ii);
        PVG_Z_bottom(ii)= Z_p(end); PVG_Vs_bottom(ii)= Vs_p(end);
        PVG_Z(ii)= (PVG_Z_top(ii)+PVG_Z_bottom(ii))/2;
        iz= find(abs(Z-PVG_Z(ii))==min(abs(Z-PVG_Z(ii)))); iz= iz(1);
        PVG_Vs(ii)= Vs(iz);
        PVG_amp(ii)= (PVG_Vs(ii) - NVG2_Vs_bottom(ii))/mean([PVG_Vs(ii) NVG2_Vs_bottom(ii)]);
    end
    
    
% Lowest positive peak is above last negative peak. 
elseif Z_p(end) < Z_np(end)
    % Sub-cases:
    % There is only one NVG
    if length(Z_np(Z_np>zmoh+10)) == 1
        NVG1_Z_top(ii)= Z_p(end); NVG1_Vs_top(ii)= Vs_p(end);
        NVG1_Z_bottom(ii)= Z_np(end); NVG1_Vs_bottom(ii)= Vs_np(end);
        
        % NVG center info
        NVG1_Z(ii)= (Z_np(end)+Z_p(end))/2;
        iz= find(abs(Z-NVG1_Z(ii))==min(abs(Z-NVG1_Z(ii)))); iz= iz(1);
        NVG1_Vs(ii)= Vs(iz);
        
        % Again, a few special cases.
        if ii==462
            % This one still has a PVG
            PVG_Z(ii)= (172+254)/2;
            iz= find(abs(Z-PVG_Z(ii))==min(abs(Z-PVG_Z(ii)))); iz= iz(1);
            PVG_Vs(ii)= Vs(iz);
            
            PVG_Z_top(ii)= 172;
            iz= find(abs(Z-PVG_Z_top(ii))==min(abs(Z-PVG_Z_top(ii)))); iz= iz(1);
            PVG_Vs_top(ii)= Vs(iz);
            PVG_Z_bottom(ii)= 254;
            iz= find(abs(Z-PVG_Z_bottom(ii))==min(abs(Z-PVG_Z_bottom(ii)))); iz= iz(1);
            PVG_Vs_bottom(ii)= Vs(iz);
            
            PVG_amp(ii)= (PVG_Vs(ii) - NVG1_Vs_bottom(ii))/mean([PVG_Vs(ii) NVG1_Vs_bottom(ii)]);
        end
    
    else % Multiple NVGs
        gVs = gradient(Vs,Z);
        
        ipt= find(Z_p<Z_np(end-1)); ipt= ipt(end);
        NVG1_Z_top(ii)= Z_p(ipt); NVG1_Vs_top(ii)= Vs_p(ipt);
        NVG1_Z_bottom(ii)= Z_np(end-1); NVG1_Vs_bottom(ii)= Vs_np(end-1);
        
        ipb= find(Z_p<Z_np(end)); 
        NVG2_Z_top(ii)= Z_p(end); NVG2_Vs_top(ii)= Vs_p(end);
        NVG2_Z_bottom(ii)= Z_np(end); NVG2_Vs_bottom(ii)= Vs_np(end);
        
        % NVG Center
        NVG1_Z(ii)= (Z_np(end-1)+Z_p(end-1))/2;
        iz= find(abs(Z-NVG1_Z(ii))==min(abs(Z-NVG1_Z(ii)))); iz= iz(1);
        NVG1_Vs(ii)= Vs(iz);
        
        NVG2_Z(ii)= (Z_np(end)+Z_p(end))/2;
        iz= find(abs(Z-NVG2_Z(ii))==min(abs(Z-NVG2_Z(ii)))); iz= iz(1);
        NVG2_Vs(ii)= Vs(iz);
        
        % PVG might be between them.
        PVG_Z_top(ii)= NVG1_Z_bottom(ii); PVG_Vs_top(ii)= NVG1_Vs_bottom(ii);
        PVG_Z_bottom(ii)= Z_p(end); PVG_Vs_bottom(ii)= Vs_p(end);
        PVG_Z(ii)= (PVG_Z_top(ii)+PVG_Z_bottom(ii))/2;
        iz= find(abs(Z-PVG_Z(ii))==min(abs(Z-PVG_Z(ii)))); iz= iz(1);
        PVG_Vs(ii)= Vs(iz);
        
    end
    
    end
    %end
end



% Reliability: if in group is1a, Reliability is high. If index is2,
% Reliability is lower.
is1a= NVG1_Vs_bottom<=4.32 & NVG1_Vs_bottom>0;
is2= NVG1_Vs_bottom>4.32;

%% What is LAB vs. MLD?

LAB_Z= zeros(len,1); MLD_Z= zeros(len,1); LAB_uncertain= zeros(len,1);
LAB_Vs= zeros(len,1); MLD_Vs= zeros(len,1);

% Case 1: Vs at base <= 4.32 km/s. First NVG is definitely the LAB.
LAB_Z(is1a)= NVG1_Z(is1a);
LAB_Vs(is1a)= NVG1_Vs(is1a);
gradual= (isGradual==1);

% Case 2: Vs at base of NVG1 is greater than 4.32 km/s. Check if there are one or two layers.
Vs_cutoff= 4.5; Vs_cutoff2= 4.55;
Z_cutoff= NVG2_Z > max(NVG1_Z_bottom(NVG1_Vs_bottom<=Vs_cutoff))+10;
shallow_is_MLD= (is2) & (NVG1_Z>0) & (NVG1_Vs_top>Vs_cutoff) & (NVG1_Vs_bottom > Vs_cutoff) & ~Z_cutoff & ~gradual;
shallow_is_LAB= (is2) & (NVG1_Z>0) & (NVG1_Vs_top>Vs_cutoff) & (NVG1_Vs_bottom <= Vs_cutoff) & ~Z_cutoff & ~gradual;

%deep_is_LAB= (is2) & ~shallow_is_LAB & (NVG2_Z>0) & (NVG2_Vs_top > Vs_cutoff2) & (NVG2_Vs_bottom <= Vs_cutoff2) & ~Z_cutoff & ~gradual;
deep_is_LAB2= (is2) & (NVG2_Z>0) & (NVG2_Vs_top > Vs_cutoff2) & (NVG2_Vs_bottom <= Vs_cutoff2) & ~Z_cutoff & ~gradual;
deep_is_LAB= deep_is_LAB2; 
shallow_is_MLD(shallow_is_LAB & deep_is_LAB2)= true;
shallow_is_LAB(shallow_is_LAB & deep_is_LAB2)= false;


% Now, populate depth and Vs at center of LAB and MLD.
LAB_Z(shallow_is_LAB)= NVG1_Z(shallow_is_LAB);
LAB_Vs(shallow_is_LAB)= NVG1_Vs(shallow_is_LAB);
LAB_Z(deep_is_LAB)= NVG2_Z(deep_is_LAB);
LAB_Vs(deep_is_LAB)= NVG2_Vs(deep_is_LAB);

MLD_Z(shallow_is_MLD)= NVG1_Z(shallow_is_MLD);
MLD_Vs(shallow_is_MLD)= NVG1_Vs(shallow_is_MLD);


LAB_amp= zeros(len,1);
LAB_amp(is1a)= 100*(NVG1_Vs_top(is1a) - NVG1_Vs_bottom(is1a))./mean([NVG1_Vs_top(is1a) NVG1_Vs_bottom(is1a)],2);
LAB_amp(shallow_is_LAB)= 100*(NVG1_Vs_top(shallow_is_LAB) - NVG1_Vs_bottom(shallow_is_LAB))./mean([NVG1_Vs_top(shallow_is_LAB) NVG1_Vs_bottom(shallow_is_LAB)],2);
LAB_amp(deep_is_LAB)= 100*(NVG2_Vs_top(deep_is_LAB) - NVG2_Vs_bottom(deep_is_LAB))./mean([NVG2_Vs_top(deep_is_LAB) NVG2_Vs_bottom(deep_is_LAB)],2);



% Populate LAB/MLD width vectors
LAB_width= zeros(len,1);
LAB_width(is1a)= NVG1_Z_bottom(is1a) - NVG1_Z_top(is1a);
LAB_width(shallow_is_LAB)= NVG1_Z_bottom(shallow_is_LAB) - NVG1_Z_top(shallow_is_LAB);
LAB_width(deep_is_LAB)= NVG2_Z_bottom(deep_is_LAB) - NVG2_Z_top(deep_is_LAB);

MLD_width= zeros(len,1);
MLD_width(shallow_is_MLD)= NVG1_Z_bottom(shallow_is_MLD) - NVG1_Z_top(shallow_is_MLD);
MLD_amp(shallow_is_MLD)= 100*(NVG1_Vs_top(shallow_is_MLD) - NVG1_Vs_bottom(shallow_is_MLD))./mean([NVG1_Vs_top(shallow_is_MLD) NVG1_Vs_bottom(shallow_is_MLD)],2);
gradual= (isGradual==1 | LAB_width>150);

% Update everything to exclude gradual points
deep_is_LAB(gradual)= false; shallow_is_LAB(gradual)= false;
shallow_is_MLD(gradual)= false;
LAB_Z(gradual)= 0; LAB_width(gradual)= 0; LAB_amp(gradual)= 0;


% Basic stats: how many NVGs?
nNVG1= sum(NVG1_Z>0); 
nLAB= sum(LAB_Z>0); nMLD= sum(MLD_Z>0);


%% PVGs
isPVG= PVG_Z>LAB_Z & LAB_Z>0 & is1a;
%Percent contrast definition.
isPVG_amp= PVG_Z>LAB_Z & LAB_Z>0;% & is1a; % gradient definition.

%% Save info
Lat= mdls.lat; Lon= mdls.lon;
save data/Gradients_all.mat Lat Lon NVG1_Z NVG1_Z_bottom NVG2_Z NVG2_Z_bottom NVG1_Vs NVG1_Vs_bottom NVG2_Vs NVG2_Vs_bottom is1a Moho_Z
save data/Gradients_all.mat PVG_Z PVG_Z_top PVG_Z_bottom PVG_Vs PVG_Vs_top PVG_Vs_bottom gradual PVG_amp isPVG_amp -append
save data/Gradients_all.mat LAB_Vs LAB_width LAB_Z LAB_amp MLD_Vs MLD_width MLD_amp MLD_Z shallow_is_LAB shallow_is_MLD deep_is_LAB -append



