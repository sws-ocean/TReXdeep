%% Simple budget integration
addpath(genpath('/Users/samst/Dropbox/UBC/TReX_deepII/'));
addpath(genpath('/Users/samst/Dropbox/UBC/TReX_paper/'));
addpath(genpath('/Users/samst/Dropbox/UBC/GSW/'));
addpath(genpath('/Users/samst/Dropbox/UBC/m_map/'));
addpath(genpath('/Users/samst/Dropbox/UBC/Misc/'));

clear

load data20230707.mat

%% How much modelled tracer lef the Gulf
mTracer=load('conc_20220604_K10');
leavePct=nansum(mTracer.C(mTracer.lat<47))/nansum(mTracer.C(:))*100;

%% Calculate hypsography
% Load bathy and pull out station locs
bathy=load('GoSL_Zhires.mat');
bathy.Z(bathy.Z==-9999)=NaN;

% Decimate bathy to speed things up
dec=5;
bathy.lon=bathy.lon(1:dec:end,1:dec:end);
bathy.lat=bathy.lat(1:dec:end,1:dec:end);
bathy.Z=bathy.Z(1:dec:end,1:dec:end)*-1;

allZ=bathy.Z(:);
allA=cdtarea(bathy.lat,bathy.lon);allA=allA(:);
Az=zeros(ceil(abs(max(bathy.Z(:)))),1);
Az_z=0:max(allZ);
for i=0:max(allZ)
    Az(i+1)=sum(allA(allZ>i));
end
  
%% Convert A(z) to A(z*)
sigma_grid1_Az_star=NaN(length(Az_z),length(CTD1));
for i=1:length(CTD1)
    sigma_grid1_Az_star(:,i)=interp1(CTD1(i).depSM,CTD1(i).sigma,Az_z);
end
mnSig1=mean(sigma_grid1_Az_star,2,'omitnan');
msk=~isnan(mnSig1);
Az_star1=interp1(mnSig1(msk),Az(msk),sigCo1);

sigma_grid2_Az_star=NaN(length(Az_z),length(CTD2));
for i=1:length(CTD2)
    sigma_grid2_Az_star(:,i)=interp1(CTD2(i).depSM,CTD2(i).sigma,Az_z);
end
mnSig2=mean(sigma_grid2_Az_star,2,'omitnan');
msk=~isnan(mnSig2);
Az_star2=interp1(mnSig2(msk)+rand(size(mnSig2(msk)))/1e10,Az(msk),sigCo2);

%% BOE
GoSLv=trapz(zStar_grid,Az_star1); % m3
GoSLm=GoSLv*1026; % kg
GoSLc=(mean(horzcat(CTD1.c_ij),1:2)/1e15); % mol/kg
GoSLmm=GoSLm*GoSLc; % total mol
BOE1=GoSLmm*mm_sf5cf3 % total g

GoSLc=(mean(horzcat(CTD2.c_ij),1:2)/1e15); % mol/kg
GoSLmm=GoSLm*GoSLc; % total mol
BOE2=GoSLmm*mm_sf5cf3 % total g

%% Budget eq 1 Holtermann et al

Azi=interp1(zStar_grid-275,Az_star1,zHat_grid);
tmpC=nanmean(horzcat(CTD1.c_ij),2);
stdC=std(horzcat(CTD1.c_ij),0,2,'omitnan');
MJun22=simps(zHat_grid(~isnan(Azi)),Azi(~isnan(Azi)).*...
    tmpC(~isnan(Azi)),1)*1026.25/1e15*mm_sf5cf3

Azi=interp1(zStar_grid-275,Az_star2,zHat_grid);
tmpC=nanmean(horzcat(CTD2.c_ij),2);
stdC=std(horzcat(CTD2.c_ij),0,2,'omitnan');
MOct22=simps(zHat_grid(~isnan(Azi)),Azi(~isnan(Azi)).*...
   tmpC(~isnan(Azi)),1)*1026.25/1e15*mm_sf5cf3

% z is a depth coordinate in m
% Az is A(z), a depth variable hypsography in m2
% cz is c(z), a depth variable concentration of a tracer in femtomols/kilogram
% rho is a density conversion factor of 1026 kg/m3 for seawater
% mm_sf5cf3 is the molar mass of the tracer compound in g/mol
% simps is an integration function
% M=simps(z,Az.*cz,1)*rho/1e15*mm_sf5cf3 in g


%% Flushing time
Tt=-(CTD2(1).mtime(1)-CTD1(1).mtime(1))/log(MOct22/MJun22)/365 % yrs

%% Integrate Ij field
Ar=cdtarea(xx,yy,'km2')*1e6;%m2 
Im1=sum(zi.*Ar,1:2,'omitnan')/1e15*mm_sf5cf3*1026
Im2=sum(zi2.*Ar,1:2,'omitnan')/1e15*mm_sf5cf3*1026

%% Initial mass
M0=MOct22*(exp((CTD2(1).mtime(1)-datenum(2021,10,25,12,00,00))...
    /Tt))
M0 =MJun22/exp(-(CTD2(1).mtime(1)-datenum(2021,10,25,12,00,00))/Tt)...
    *exp(MOct22/MJun22*Tt);
M0 = MJun22 / exp(-(CTD1(1).mtime(1)-datenum(2021,10,25,12,00,00))/Tt)
M0 = MOct22 / exp(-(CTD2(1).mtime(1)-datenum(2021,10,25,12,00,00))/Tt)

%% Budget per area
lat_lim=[46 53];
lon_lim=[-67.5 -55];
m_proj('equidistant','lon',lon_lim,'lat',lat_lim);
f1=figure('units','centimeters','outerposition',...
    [0.01 0.01 20 20],'color','w');
ax1=axes; hold on;%('position',[.0 .0 .1 1]); hold on;
m_proj('equidistant','lon',lon_lim,'lat',lat_lim);
[c,h]=m_contour(bathy.lon,bathy.lat,bathy.Z,0:-100:-500,'color',...
    [0.95 0.95 0.95],'linewidth',0.5);
clabel(c,h,'LabelSpacing',500,'color',[0.95 0.95 0.95],'fontsize',6);
m_gshhs_i('patch',[0.9 0.9 0.9],'edgecolor','none');
m_grid('linestyle','none','linewidth',0.1,'tickdir','out',...
    'xaxisloc','top','yaxisloc','right','fontsize',6);
s=m_scatter(dLon,dLat,200,'p','filled',...
    'markerfacecolor',[0.929,0.694,0.125],'markeredgecolor','k');
m_contour(Zlon,Zlat,elev,[-275 -275],'color',[0.2 0.2 0.2]);
m_scatter(stlon1,stlat1,25,'filled',...
    'markerfacecolor','b','markeredgecolor','w','markerfacealpha',0.3);
m_scatter(stlon2,stlat2,25,'filled',...
    'markerfacecolor','r','markeredgecolor','w','markerfacealpha',0.3);

LC.lon=[-59.2654  -61.6800  -63.2289  -57.8531 -56.5091  -58.0125 ...
     -60.6321  -59.6526  -59.0831  -58.4453  -57.9214];
LC.lat=[48.5151   49.0887   49.8434   51.6999 51.1112   49.5415 ...
    49.5319   49.2753   49.8036   50.1961   50.7545];

LC.station_bound=alphaShape(LC.lon',LC.lat',1,'HoleThreshold',15);
[l,ll]=m_ll2xy(LC.lon,LC.lat);
tmpshp=alphaShape(l',ll',1,'HoleThreshold',15);
LC.AS=plot(tmpshp,'facealpha',0.1,'edgecolor','none','facecolor',rgb_x('forest green'));

idx1=find(inShape(LC.station_bound,stlon1,stlat1) & [CTD1.x_dist]'>0);
idx2=find(inShape(LC.station_bound,stlon2,stlat2) & [CTD2.x_dist]'>0);

idx4=find(inShape(LC.station_bound,bathy.lon,bathy.lat));
allZ=bathy.Z(idx4);
allA=cdtarea(bathy.lat,bathy.lon);allA=allA(idx4);
Az=zeros(ceil(abs(max(bathy.Z(:)))),1);
Az_z=0:ceil(max(abs(bathy.Z(:))))-1;
for i=0:max(allZ)
    Az(i+1)=sum(allA(allZ>i));
end

tmp1=horzcat(CTD1.c_ij);
tmp2=horzcat(CTD2.c_ij);

Azi=interp1(Az_z'-275,Az,zHat_grid);
MLC1=simps(zHat_grid,Azi.*mean(tmp1(:,idx1),2),1)*1026.25/1e15*mm_sf5cf3
MLC2=simps(zHat_grid,Azi.*mean(tmp2(:,idx2),2),1)*1026.25/1e15*mm_sf5cf3


% sanity check
idx1=find(~inShape(LC.station_bound,stlon1,stlat1) & [CTD1.x_dist]'>0);
idx2=find(~inShape(LC.station_bound,stlon2,stlat2) & [CTD2.x_dist]'>0);

idx4=find(~inShape(LC.station_bound,bathy.lon,bathy.lat));
allZ=bathy.Z(idx4);
allA=cdtarea(bathy.lat,bathy.lon);allA=allA(idx4);
Az=zeros(ceil(abs(max(bathy.Z(:)))),1);
Az_z=0:ceil(max(abs(bathy.Z(:))))-1;
for i=0:max(allZ)
    Az(i+1)=sum(allA(allZ>i));
end

tmp1=horzcat(CTD1.c_ij);
tmp2=horzcat(CTD2.c_ij);

Azi=interp1(Az_z'-275,Az,zHat_grid);
MLC1=simps(zHat_grid,Azi.*mean(tmp1(:,idx1),1:2),1)*1026.25/1e15*mm_sf5cf3
MLC2=simps(zHat_grid,Azi.*mean(tmp2(:,idx2),1:2),1)*1026.25/1e15*mm_sf5cf3



%% Flushing time
ff=fit([1;CTD1(1).mtime(1)-datenum(2021,10,25,12,00,00);...
    CTD2(1).mtime(1)-datenum(2021,10,25,12,00,00); 600],[5000;MJun22;MOct22;50],...
    'exp1');
ft=feval(ff,1:365*5);
coeffvals=coeffvalues(ff);
T_f = -log(MOct22/MJun22)/coeffvals(2);




% 1. Calculate the ratio of the tracer mass in October to the tracer mass in July. Let's call this ratio R.
% 
% 2. Assume that the tracer concentration in the estuary decays exponentially over time, with a decay rate constant k. Then, the ratio R can be related to the flushing time T_f of the estuary by the following equation:
% 
% R = exp(-k * T_f)
% 
% 3. Rearrange the equation to solve for the flushing time T_f:
% 
% T_f = -ln(R) / k
% 
% In other words, the flushing time of the estuary is given by the natural 
% logarithm of the ratio R, divided by the decay rate constant k. 
% The decay rate constant k can be estimated from the tracer measurements 
% by fitting an exponential decay curve to the data. 

%% Budget eq 1 Holtermann et al but only for tainted region

m_proj('equidistant','lon',lon_lim,'lat',lat_lim);
zz1=interp2(xx,yy,zi,bathy.lon,bathy.lat);
idx=find(zz1>2 & bathy.Z>150);
b=boundary(bathy.lon(idx),bathy.lat(idx),1);
as=alphaShape(bathy.lon(idx),bathy.lat(idx),1,'HoleThreshold',15);

allZ=bathy.Z(idx);
Az=zeros(ceil(abs(max(allZ))),1);
Az_z=0:max(allZ);
for i=0:max(allZ)
    Az(i+1)=sum(allA(allZ>i));
end

IS=inShape(as,stlon1,stlat1);
tmp=horzcat(CTD1.c_ij);

Azi=interp1(Az_z'-275,Az,zHat_grid);
idx=~isnan(Azi);
GoSLci=(mean(tmp(:,IS),2)/1e15); % mol/kg
MJun22=simps(zHat_grid(idx),Azi(idx).*GoSLci(idx),1)*1026.25/1e15*mm_sf5cf3
%%

JunArea=sum(allA(idx))

figure; hold on
m_contour(Zlon,Zlat,elev,[-25 -50 -75],'color',[0.7 0.7 0.7]);
[c,h]=m_contourf(bathy.lon,bathy.lat,zz1,0:1:40,'linestyle','none');
m_scatter(bathy.lon(idx),bathy.lat(idx),20);
% plot(tmpshp,'facealpha',0.1,'facecolor','r','edgecolor','none');
m_plot(bathy.lon(idx(b)),bathy.lat(idx(b)),'k','linewidth',2);
caxis([0 40]);
% m_scatter(bathy.lon,bathy.lat,20);
colormap(gca,seminfhaxby);
m_gshhs_i('patch',rgb('light grey'),'edgecolor','none');
m_grid('linestyle','none','linewidth',0.1,'tickdir','out',...
'xaxisloc','bottom','yaxisloc','left','fontsize',12);

zz2=interp2(xx,yy,zi2,bathy.lon,bathy.lat);
idx=find(zz2>2 & bathy.Z>200);
b=boundary(bathy.lon(idx),bathy.lat(idx),1);
OctArea=sum(allA(idx))

figure; hold on
m_contour(Zlon,Zlat,elev,[-25 -50 -75],'color',[0.7 0.7 0.7]);
[c,h]=m_contourf(bathy.lon,bathy.lat,zz2,0:1:40,'linestyle','none');
m_scatter(bathy.lon(idx),bathy.lat(idx),20);
% plot(tmpshp,'facealpha',0.1,'facecolor','r','edgecolor','none');
m_plot(bathy.lon(idx(b)),bathy.lat(idx(b)),'k','linewidth',2);
caxis([0 40]);
% m_scatter(bathy.lon,bathy.lat,20);
colormap(gca,seminfhaxby);
m_gshhs_i('patch',rgb('light grey'),'edgecolor','none');
m_grid('linestyle','none','linewidth',0.1,'tickdir','out',...
'xaxisloc','bottom','yaxisloc','left','fontsize',12);


% as=alphaShape(bathy.lon(zz>0.5),bathy.lat(zz>0.5),1,'HoleThreshold',15);
% [l,ll]=m_ll2xy(bathy.lon,bathy.lat);
% tmpshp=alphaShape(l(zz>0.5),ll(zz>0.5),1,'HoleThreshold',15);
%%
Hspread=(OctArea-JunArea)/((CTD2(1).mtime(1)-CTD1(1).mtime(1))*24*60*60)
