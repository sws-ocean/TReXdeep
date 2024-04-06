%% Results figure from other plotting programs
addpath(genpath('/Users/samst/Dropbox/UBC/TReX_paper/'));
addpath(genpath('/Users/samst/Dropbox/UBC/TReX_deepII/'));
addpath(genpath('/Users/samst/Dropbox/UBC/GSW/'));
addpath(genpath('/Users/samst/Dropbox/UBC/m_map/'));
addpath(genpath('/Users/samst/Dropbox/UBC/Misc/'));
addpath(genpath('/Users/samst/Dropbox/UBC/drifter/'));

clear

load data20230707.mat

model10=load('model10_z.mat');

% Create model fields
model10.Cp=squeeze(sum(sum(model10.CDzi,'omitnan')));
model10.Cz=squeeze(sum(sum(model10.CDzi,'omitnan')));
model10.z=[1:500]';

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
% zHat_grid=interp1(mnSig1(msk),Az(msk),sigCo1);
Az_star1=interp1(mnSig1(msk)+rand(size(mnSig1(msk)))/1e10,Az(msk),sigCo1);
[~,idx1]=min(abs(mnSig1-27.26));
Az_hat1=interp1(zStar_grid-zStar_grid(idx1),Az_star1,zHat_grid);

sigma_grid2_Az_star=NaN(length(Az_z),length(CTD2));
for i=1:length(CTD2)
    sigma_grid2_Az_star(:,i)=interp1(CTD2(i).depSM,CTD2(i).sigma,Az_z);
end
mnSig2=mean(sigma_grid2_Az_star,2,'omitnan');
msk=~isnan(mnSig2);
Az_star2=interp1(mnSig2(msk)+rand(size(mnSig2(msk)))/1e10,Az(msk),sigCo2);
[~,idx2]=min(abs(mnSig2-27.26));
Az_hat2=interp1(zStar_grid-zStar_grid(idx2),Az_star2,zHat_grid);

%% Load in float info 
T1=load('TReXIfloat.mat');
T2=load('TReXIIfloat.mat');

% Average popup position and velocity
for i=1:length(T1.F)
    T1.F(i).bearing=azimuth(T1.F(i).deploymentLat,T1.F(i).deploymentLon,...
        T1.F(i).resurfLat,T1.F(i).resurfLon);
end
bearings=[T1.F(i).bearing];
mean_sin = mean(sin(deg2rad(bearings)));
mean_cos = mean(cos(deg2rad(bearings)));
avg_bearing_deg = rad2deg(atan2(mean_sin, mean_cos));
T1_deg = mod(avg_bearing_deg, 360);
T1v=mean([T1.F.speedCMsec]); %cm s-1
T1s=std([T1.F.speedCMsec]);
T1u=T1s/sqrt(length([T1.F.speedCMsec]));
distance_deg =T1v*0.036/111.11*24*30;
delta_lat = distance_deg * cos(deg2rad(avg_bearing_deg));
delta_long = distance_deg * sin(deg2rad(avg_bearing_deg)) /...
    cos(deg2rad(T1.F(i).deploymentLat));
latT1=T1.F(i).deploymentLat + delta_lat;
lonT1=T1.F(i).deploymentLon + delta_long;

for i=1:length(T2.F)
    T2.F(i).bearing=azimuth(T2.F(i).deploymentLat,T2.F(i).deploymentLon,...
        T2.F(i).resurfLat,T2.F(i).resurfLon);
end
bearings=[T2.F(i).bearing];
mean_sin = mean(sin(deg2rad(bearings)));
mean_cos = mean(cos(deg2rad(bearings)));
avg_bearing_deg = rad2deg(atan2(mean_sin, mean_cos));
T2_deg = mod(avg_bearing_deg, 360);
T2v=mean([T2.F.speedCMsec]); %cm s-1
T2s=std([T2.F.speedCMsec]);
T2u=T2s/sqrt(length([T2.F.speedCMsec]));
distance_deg =T2v*0.036/111.11*24*30;
delta_lat = distance_deg * cos(deg2rad(avg_bearing_deg));
delta_long = distance_deg * sin(deg2rad(avg_bearing_deg)) /...
    cos(deg2rad(T2.F(i).deploymentLat));
latT2=T2.F(i).deploymentLat + delta_lat;
lonT2=T2.F(i).deploymentLon + delta_long;

% diffusivity? 
% T1
k=nchoosek(1:length(T1.F),2);
count=0; l=[];
for i=1:length(k)
    if abs(T1.F(k(i,1)).dTday-T1.F(k(i,2)).dTday)<3
        count=count+1;
        l(count,1:2)=[gsw_distance([...
            T1.F(k(i,1)).resurfLon T1.F(k(i,2)).resurfLon],...
            [T1.F(k(i,1)).resurfLat T1.F(k(i,2)).resurfLat]) ...
            mean([T1.F(k(i,1:2)).dTday])];
    end
end

% T2
k=nchoosek(1:length(T2.F),2);
for i=1:length(k)
    if abs(T2.F(k(i,1)).dTday-T2.F(k(i,2)).dTday)<3
        count=count+1;
        l(count,1:2)=[gsw_distance([...
            T2.F(k(i,1)).resurfLon T2.F(k(i,2)).resurfLon],...
            [T2.F(k(i,1)).resurfLat T2.F(k(i,2)).resurfLat]) ...
            mean([T2.F(k(i,1:2)).dTday])];
    end
end

% pair separation velocity m/s
l(:,3)=l(:,1)./(l(:,2)*86400);

% add 0 disp at time 0? 
doubleK=0.5*mean(l(:,1).^2./(l(:,2)*86400))

I1=bootrnd(size(l(:,1),1),1000);
doubleKbt=[];
for j=1:size(I1,2)
    doubleKbt(j,1)=0.5*mean(l(I1(:,j),1).^2./(l(I1(:,j),2)*86400));
end

% std(doubleKbt)/sqrt(




%% Plotting
f1=figure('units','centimeters','outerposition',...
    [0.01 0.01 19 15],'color','w');
col=[255 214 140]/255; % YELLOW!

lat_lim=[46.75 51.5];
lon_lim=[-69.5 -57.5];

% Horizontal %
lev=0:1e3:30e3;
cax=[0 30e3];

%%%%%% maps %%%%%%
ax1=axes('position',[0.1 0.5 0.45 0.45]); hold on;
m_proj('equidistant','lon',lon_lim,'lat',lat_lim);
m_contour(Zlon,Zlat,elev,[-25 -50 -75],'color',[0.7 0.7 0.7]);
[c,h]=m_contourf(xx,yy,zi.*1026,lev,'linestyle','none');% This includes a density conversion factor to convert from fmol/kg/m to fmol/m2
caxis(cax);
tmp=seminfhaxby;tmp2=linspecer;%.*repmat(1-.2*cos([0:127]'/64*2*pi).^4,1,3);
colormap(gca,flipud([flipud(tmp2(1:end,:));tmp(1,:)]));
m_gshhs_f('patch',col,'edgecolor',[0.3 0.3 0.3]);
dLat=48.244; dLon=-59.952;

s=m_scatter(dLon,dLat,100,'p','filled',...
    'markerfacecolor',rgb_x('light violet'),'markeredgecolor','k',...
    'linewidth',0.3);

m=max([T1.F.dTday]);
for i=1:length(T1.F)
    s=m_scatter(T1.F(i).resurfLon,T1.F(i).resurfLat,50-(T1.F(i).dTday/m*30)...
        ,'filled','markerfacecolor',rgb_x('off white'),'markeredgecolor','k',...
        'linewidth',0.4,'markerfacealpha',0.8);
end

m_scatter(-62,51,50-T1.F(8).dTday/m*30,...
    'filled','markerfacecolor',rgb_x('off white'),'markeredgecolor','k',...
        'linewidth',0.4,'markerfacealpha',0.8)
m_scatter(-62,50.8,50-T1.F(1).dTday/m*30,...
    'filled','markerfacecolor',rgb_x('off white'),'markeredgecolor','k',...
        'linewidth',0.4,'markerfacealpha',0.8)
m_text(-61.75,51.2,'Float popups','fontsize',6,'fontweight','bold');
m_text(-61.75,51,'1 week','fontsize',6);
m_text(-61.75,50.8,'4 weeks','fontsize',6);

% cm=cmocean('gray');
% weig=round(length(cm).*([T1.F.dTday]/max([T1.F.dTday])));
% for i=1:length(T1.F)
%     s(i)=m_scatter(T1.F(i).resurfLon,T1.F(i).resurfLat,50,...
%         'filled','markerfacecolor',cm(weig(i),:),'markeredgecolor',...
%         'none');
% end

m_scatter(stlon1,stlat1,5,'x','markeredgecolor',[.5 .5 .5],...
    'linewidth',0.5);
m_grid('linestyle','none','linewidth',0.1,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','left','fontsize',12,'xticks',[]);
set(gca,'clipping','on');
m_text(-69.3,51,'a)','fontsize',8);
m_plot([-61.6730, -59.2094],[49.1034,48.4974],'--','color',[.4 .4 .4]);

hh=m_contfbar(gca,0.23,[0.125 0.4],c,h,'endpiece','yes','xaxisloc','top',...
    'axfrac',0.04,'edgecolor','none');
ylabel(hh,{'\it{I_j} \rm(fmol m^{-2})'});


%%%%%%%%%%%%%%%
ax2=axes('position',[0.1 0.05 0.45 0.45]); hold on;
m_proj('equidistant','lon',lon_lim,'lat',lat_lim);
m_contour(Zlon,Zlat,elev,[-25 -50 -75],'color',[0.7 0.7 0.7]);
[c,h]=m_contourf(xx,yy,zi2.*1026,lev,'linestyle','none');% This includes a density conversion factor to convert from fmol/kg/m to fmol/m2
clabel(c,h,'LabelSpacing',500,'color',[0.9 0.9 0.9],'fontsize',6);

tmp=seminfhaxby;tmp2=linspecer;%.*repmat(1-.2*cos([0:127]'/64*2*pi).^4,1,3);
colormap(gca,flipud([flipud(tmp2(1:end,:));tmp(1,:)]));
caxis(cax); 
m_gshhs_f('patch',col,'edgecolor',[0.3 0.3 0.3]);
s=m_scatter(dLon,dLat,100,'p','filled',...
    'markerfacecolor',rgb_x('light violet'),'markeredgecolor','k',...
    'linewidth',0.3);
s=m_scatter(T2.F(1).deploymentLon,T2.F(1).deploymentLat,30,...
    'x','markeredgecolor','k','linewidth',2);

m=max([T2.F.dTday]);
for i=1:length(T2.F)
    s=m_scatter(T2.F(i).resurfLon,T2.F(i).resurfLat,50-(T2.F(i).dTday/m*30)...
        ,'filled','markerfacecolor',rgb_x('off white'),'markeredgecolor','k',...
        'linewidth',0.4,'markerfacealpha',0.8);
end
m_scatter(stlon2,stlat2,5,'x','markeredgecolor',[.5 .5 .5],...
    'linewidth',0.5);

m_grid('linestyle','none','linewidth',0.1,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','left','fontsize',12);
set(gca,'clipping','on');
m_text(-69.3,51,'d)','fontsize',8);
m_plot([-61.6730, -59.2094],[49.1034,48.4974],'--','color',[0.4 .4 .4]);
m_ruler([0.05 0.225],0.125,'linewidth',0.5);

m_annotation('textarrow',[dLon-1.25 dLon-0.1],[dLat-0.4 dLat],...
    'string',{'Injection';'site'},'fontsize',6,'Fontweight','normal',...
    'horizontalalignment','center','headlength',5);

% Vertical %
ax3=axes('position',[0.625 0.525 0.15 0.425]); hold on;
% ll=plot((model10.Cp/max(model10.Cp))*0.8,model10.zco+292,'LineWidth',2,...
%     'Color',rgb_x('light blue'));
% ll.Color=[ll.Color 0.2];
% fill(maxC1zi,maxZ,[0.9 0.9 0.9],'linestyle','none');
s1=scatter(vertcat(CTD1.TR_SF5),vertcat(CTD1.zStar),10,'filled',...
    'markerfacecolor',rgb_x('blue'),'markeredgecolor','none',...
    'markerfacealpha',0.3');

s2=scatter(600,600,10,'filled',...
    'markerfacecolor',rgb_x('red'),'markeredgecolor','none',...
    'markerfacealpha',0.4');

axis ij; grid on;
ylabel('\it{z^*} \rm{(m)}');
xlim([0 0.9]); ylim([0 500]);
xticks(0:.25:1);
ax3.XTickLabel=[];
text(0.1,20,'b)','fontsize',8);
line([0 0.9],[zStar_grid(sigIdx1) zStar_grid(sigIdx1)],...
    'color','k','linestyle',':','linewidth',1);l2=plot(NaN,NaN,'color',[0.4 0.4 0.4]);
l3=plot(NaN,NaN,'--','color',[0.4 0.4 0.4]);

ax3i=axes('position',[0.625 0.525 0.15 0.425],'color','none'); hold on;
plot(Az_star1,zStar_grid,'color',[0.4 0.4 0.4]);
plot(Az_star1(1)*exp(-zStar_grid/200),zStar_grid,'--','color',[0.4 0.4 0.4]);
ylim([0 500]);yticks([]); xlim([0 max([Az_star1;Az_star2])])
xticklabels({'0','1×10^{11}','2×10^{11}'});xtickangle(0);
set(gca,'xaxislocation','top','ydir','reverse',...
    'xcolor',[0.4 0.4 0.4],'TickLength',[0.03, 0.01]);

% Create leg 2 datasets
tmpMsk=logical(vertcat(CTD2.isMir));
tmpC=vertcat(CTD2.TR_SF5Mir);tmpCMir=tmpC(tmpMsk);tmpCnotMir=tmpC(~tmpMsk);
tmpZ=vertcat(CTD2.zStarMir);tmpZMir=tmpZ(tmpMsk);tmpZnotMir=tmpZ(~tmpMsk);
tmpS=vertcat(CTD2.sigMir);tmpSMir=tmpS(tmpMsk);tmpSnotMir=tmpS(~tmpMsk);

ax4=axes('position',[0.8 0.525 0.15 0.425]); hold on;
% fill(maxC2zi,maxZ,[0.9 0.9 0.9],'linestyle','none');
scatter(tmpCnotMir,tmpZnotMir,10,'filled',...
    'markerfacecolor',rgb_x('red'),'markeredgecolor','none',...
    'markerfacealpha',0.4');
axis ij; grid on;
xlim([0 0.9]); ylim([0 500]);
xticks(0:.25:1);
ax4.YTickLabel=[];
ax4.XTickLabel=[];
text(0.1,20,'c)','fontsize',8);
l1=line([0 0.9],[zStar_grid(sigIdx2) zStar_grid(sigIdx2)],...
    'color','k','linestyle',':','linewidth',1);
% l.Color=[l.Color 0.7];ll.Color=[ll.Color 0.7];
l2=plot(NaN,NaN,'color',[0.4 0.4 0.4]);
l3=plot(NaN,NaN,'--','color',[0.4 0.4 0.4]);

ax4i=axes('position',[0.8 0.525 0.15 0.425],'color','none'); hold on;
plot(Az_star2,zStar_grid,'color',[0.4 0.4 0.4]);
plot(Az_star2(1)*exp(-zStar_grid/200),zStar_grid,'--','color',[0.4 0.4 0.4]);
ylim([0 500]);yticks([]); xlim([0 max([Az_star1;Az_star2])])
xticklabels({'0','1×10^{11}','2×10^{11}'});xtickangle(0);
set(gca,'xaxislocation','top','ydir','reverse',...
    'xcolor',[0.4 0.4 0.4],'TickLength',[0.03, 0.01]);
xlabel('\it{A(z^*)} \rm{[m^2]}');

ax5=axes('position',[0.625 0.075 0.15 0.425]); hold on;
% fill(maxC1si,maxSig,[0.9 0.9 0.9],'linestyle','none');
scatter(vertcat(CTD1.TR_SF5),vertcat(CTD1.TR_sigma),10,'filled',...
    'markerfacecolor',rgb_x('blue'),'markeredgecolor','none',...
    'markerfacealpha',0.3');
axis ij; grid on;
ylabel('\sigma_\theta (kg m^{-3})');
xlim([0 0.9]); ylim([26.6 27.6]);
xticks(0:.25:1);
text(0.1,26.65,'e)','fontsize',8);
line([0 0.9],[27.26 27.26],'color','k','linestyle',':','linewidth',1);

ax6=axes('position',[0.8 0.075 0.15 0.425]); hold on;
% fill(maxC2si,maxSig,[0.9 0.9 0.9],'linestyle','none');
scatter(tmpCnotMir,tmpSnotMir,10,'filled',...
    'markerfacecolor',rgb_x('red'),'markeredgecolor','none',...
    'markerfacealpha',0.4');
axis ij; grid on;
xlim([0 0.9]); ylim([26.6 27.6]);
xticks(0:.25:1);
xlabel('SF_5CF_3 (fmol kg^{-1})');
ax6.YTickLabel=[];
text(0.1,26.65,'f)','fontsize',8);
line([0 0.9],[27.26 27.26],'color','k','linestyle',':','linewidth',1);

set(findall(gcf,'-property','fontsize'),'fontsize',8);
axes(ax4);
legend([l1;l2;l3],'Injection','\it{A(z^*)}',...
    '\it{A_0e^{-z^{*}/L}}','location','southeast','fontsize',6);
set(findall(gcf,'-property','fontname'),'fontname','Latin Modern Roman');
uistack(ax4,'bottom');

%% 
export_fig /Users/samst/Dropbox/UBC/TReX_paper/figures/R1/resultsR1.pdf -dpdf -nofontswap