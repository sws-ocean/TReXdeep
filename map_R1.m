%% Results figure from other plotting programs
addpath(genpath('/Users/samst/Dropbox/UBC/TReX_paper/'));
addpath(genpath('/Users/samst/Dropbox/UBC/TReX_deepII/'));
addpath(genpath('/Users/samst/Dropbox/UBC/GSW/'));
addpath(genpath('/Users/samst/Dropbox/UBC/m_map/'));
addpath(genpath('/Users/samst/Dropbox/UBC/Misc/'));

clear

%% CTD DATA
load('data20230707.mat', 'CTD1','CTD2','stBath1')
mTracer=load('conc_20220604_K10');

% Load bathy and pull out station locs
bathy=load('GoSL_Zhires.mat');
bathy.Z(bathy.Z==-9999)=NaN;
% Decimate bathy to speed things up
dec=20;
bathy.lon=bathy.lon(1:dec:end,1:dec:end);
bathy.lat=bathy.lat(1:dec:end,1:dec:end);
bathy.Z=bathy.Z(1:dec:end,1:dec:end);

stlon1=[];stlat1=[];
count=0;
for i=1:length(CTD1)
    if ~isempty(CTD1(i).TR_SF5)
        count=count+1;
        stlon1(count)=CTD1(i).longitude(1);
        stlat1(count)=CTD1(i).latitude(1);
    end
end

stlon2=[];stlat2=[];
count=0;
for i=1:length(CTD2)
    if ~isempty(CTD2(i).TR_SF5)
        count=count+1;
        stlon2(count)=CTD2(i).longitude(1);
        stlat2(count)=CTD2(i).latitude(1);
        stBath2(count)=griddata(bathy.lon,bathy.lat,bathy.Z,...
            stlon2(count),stlat2(count));
    end
end

% Load bathy and pull out station locs
bathy=load('GoSL_Zhires.mat');
bathy.Z(bathy.Z==-9999)=NaN;
% Decimate bathy to speed things up
dec=4;
bathy.lon=bathy.lon(1:dec:end,1:dec:end);
bathy.lat=bathy.lat(1:dec:end,1:dec:end);
bathy.Z=bathy.Z(1:dec:end,1:dec:end);


%% Climatological CTD data
% all_ctd=wod_rd('ocldb1625120301.1382.CTD','params',1:3);
all_ctd=wod_rd('GoSL_clim_2015_2022.CTD','params',1:3);
% tmp=wod_rd('GoSL_DOclim_2005_2022.CTD','params',3);
a=sw_dens(all_ctd.sal,all_ctd.temp,all_ctd.depth);
b=gsw_p_from_z(-all_ctd.depth,all_ctd.latitude);
SA=gsw_SA_from_SP(all_ctd.sal,b,all_ctd.longitude,all_ctd.latitude);
SR=gsw_SR_from_SP(all_ctd.sal);
CT=gsw_CT_from_t(SA,all_ctd.temp,b);
sigma=gsw_sigma0(SA,CT);

GoSL_rho=[];all_ctd.pr=[];GoSL_temp=[];GoSL_sigma=[];GoSL_SA=[];GoSL_DO=[];
GoSL_CT=[];

for i=1:length(all_ctd.mtime)
    try
        msk=~isnan(all_ctd.depth(:,i));
        GoSL_rho(:,i)=interp1(all_ctd.depth(msk,i),a(msk,i),1:400);
        GoSL_sigma(:,i)=interp1(all_ctd.depth(msk,i),sigma(msk,i),1:400);
        GoSL_temp(:,i)=interp1(all_ctd.depth(msk,i),all_ctd.temp(msk,i),1:400);
        GoSL_SR(:,i)=interp1(all_ctd.depth(msk,i),SR(msk,i),1:400);
        GoSL_CT(:,i)=interp1(all_ctd.depth(msk,i),CT(msk,i),1:400);
        all_ctd.pr(:,i)=interp1(all_ctd.depth(msk,i),b(msk,i),1:400);
    catch
        GoSL_rho(:,i)=NaN(1,400);
        GoSL_sigma(:,i)=NaN(1,400);
        GoSL_temp(:,i)=NaN(1,400);
        GoSL_SA(:,i)=NaN(1,400);
        GoSL_CT(:,i)=NaN(1,400);
        all_ctd.pr(:,i)=NaN(1,400);
    end
end

GoSL_month=month(all_ctd.mtime); GoSL_year=year(all_ctd.mtime);
GoSLNovMsk=GoSL_month>=10 & GoSL_month<=12;

GoSLmnSigma=nanmean(GoSL_sigma(:,GoSLNovMsk),2); GoSLstdSigma=std(GoSL_sigma(:,GoSLNovMsk),0,2,'omitnan');
GoSLmnTemp=nanmean(GoSL_temp(:,GoSLNovMsk),2);GoSLstdTemp=std(GoSL_temp(:,GoSLNovMsk),0,2,'omitnan');
GoSLmnSR=nanmean(GoSL_SR(:,GoSLNovMsk),2); GoSLstdSR=std(GoSL_SR(:,GoSLNovMsk),0,2,'omitnan');
% GoSLmnDO=nanmean(GoSL_DO(:,GoSLNovMsk),2); GoSLstdDO=std(GoSL_DO(:,GoSLNovMsk),0,2,'omitnan');


% load GoSL_DOclim_2005_2022.mat
DOrd;% data from https://cioosatlantic.ca/erddap/info/bio_atlantic_zone_monitoring_program_ctd/index.html
DO(1:3,:)=[];DO=table2struct(DO);

for i=1:length(DO)
    str=char(DO(i).time);
    DO(i).mtime=datenum(str2num(str(1:4)),str2num(str(6:7)),str2num(str(9:10)),...
        str2num(str(12:13)),str2num(str(15:16)),str2num(str(18:19)));
end
DOmon=month(vertcat(DO.mtime));
DONovMsk=DOmon>=10 & DOmon<=12;
DO(~DONovMsk)=[];

z=vertcat(DO.depth);
o=vertcat(DO.DOXYZZ01)*44.661;
t=vertcat(DO.TEMPPR01);
s=gsw_SR_from_SP(vertcat(DO.PSLTZZ01));

GoSLmnDO=NaN(400,1);GoSLstdDO=NaN(400,1);
GoSLmnTemp=NaN(400,1);GoSLstdTemp=NaN(400,1);
GoSLmnSR=NaN(400,1);GoSLstdSR=NaN(400,1);

for i=1:400
    GoSLmnDO(i,1)=nanmean(o(z>i & z<i+1));
    GoSLstdDO(i,1)=nanstd(o(z>i & z<i+1));
    
    GoSLmnTemp(i,1)=nanmean(t(z>i & z<i+1));
    GoSLstdTemp(i,1)=nanstd(t(z>i & z<i+1));
    
    GoSLmnSR(i,1)=nanmean(s(z>i & z<i+1));
    GoSLstdSR(i,1)=nanstd(s(z>i & z<i+1));
end


%% OTIS DATA

load OTISctd.mat

in_msk=find(inject.Pump1_SF6==1 & inject.Longitude<-50 |...
    inject.Pump2_SF6==1 & inject.Longitude<-50);
mnSig=mean(inject.Sigma0_2(in_msk));

% Find nearby CTDs? 
load TReX_injectCTD.mat

stlon=NaN(length(CTD),1);stlat=NaN(length(CTD),1);
for i=1:length(CTD)
    stlon(i)=CTD(i).longitude(1);
    stlat(i)=CTD(i).latitude(1);  
    CTD(i).sr=gsw_SR_from_SP(CTD(i).sal00);
end

%% Mapping
lat_lim=[46.75 51.5];
lon_lim=[-69.75 -57.5];
col=[255 214 140]/255; % YELLOW!
    
f1=figure('units','centimeters','outerposition',...
    [0.01 0.01 19 10],'color','w');
ax1=axes('position',[0.1 0.1 0.5 0.9]); hold on;

%%%%%% maps %%%%%%
hold on;
m_proj('equidistant','lon',lon_lim,'lat',lat_lim);
r=m_rectangle(-69.4068,47.6526,2,2,0,...
    'edgecolor','none','facecolor',rgb_x('forest green'));r.FaceColor=[r.FaceColor 0.1];
m_gshhs_f('patch',col,'edgecolor',[0.3 0.3 0.3],'linewidth',0.3);
m_grid('linestyle','none','linewidth',0.1,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','left','fontsize',8);

dLat=48.244; dLon=-59.952;
m_contour(bathy.lon,bathy.lat,bathy.Z,[-50 -100 -150 -200 -300 -400],'color',[0.9 0.9 0.9]);
% [c,h]=m_contourf(mTracer.lon,mTracer.lat,mTracer.C,0:.1e-4:2.5e-4,...
%     'linestyle','none');
tmp=seminfhaxby;tmp2=linspecer;%.*repmat(1-.2*cos([0:127]'/64*2*pi).^4,1,3);
colormap(gca,flipud([flipud(tmp2(1:end,:));tmp(1,:)]));
% caxis([0 2.75e-4]);
m_contour(bathy.lon,bathy.lat,bathy.Z,[-275 -275],...
    'color',rgb_x('dark grey'),'linewidth',0.66,'linestyle','-');
% hh=m_contfbar(gca,0.175,[0.135 0.4],c,h,'endpiece','yes','xaxisloc','top',...
%     'axfrac',0.04,'edgecolor','none');

% ylabel(hh,{'Relative';'tracer conc.'});
s=m_scatter(dLon,dLat,200,'p','filled',...
    'markerfacecolor',rgb_x('light violet'),'markeredgecolor','k');

m_scatter(stlon1(stBath1>-275),stlat1(stBath1>-275),15,'o','filled',...
    'markeredgecolor','none','markerfacecolor',rgb_x('cornflower'),'markerfacealpha',0.8);
m_scatter(stlon1(stBath1<-275),stlat1(stBath1<-275),15,'^','filled',...
    'markeredgecolor','none','markerfacecolor',rgb_x('cornflower'),'markerfacealpha',0.8);

m_scatter(stlon2(stBath2>-275),stlat2(stBath2>-275),15,'o','filled',...
    'markeredgecolor','none','markerfacecolor',rgb_x('light red'),'markerfacealpha',0.8);
m_scatter(stlon2(stBath2<-275),stlat2(stBath2<-275),15,'^','filled',...
    'markeredgecolor','none','markerfacecolor',rgb_x('light red'),'markerfacealpha',0.8);

m_annotation('textarrow',[dLon-1.25 dLon-0.1],[dLat-0.4 dLat],...
    'string',{'Injection';'site'},'fontsize',6,'Fontweight','normal',...
    'horizontalalignment','center','headlength',5);


% m_scatter(-61,51,20,'s','markerfacecolor','b','markerfacealpha',0.3,...
%     'markeredgecolor','none');
% m_scatter(-61,50.5,20,'s','markerfacecolor','r','markerfacealpha',0.3,...
%     'markeredgecolor','none');
m_text(-60.75,51,'Leg 1','fontsize',6,'color',rgb_x('cornflower'));
m_text(-60.75,50.8,'Leg 2','fontsize',6,'color',rgb_x('deep pink'));

m_text(-58.0872-0.5,48.4215-0.4,'\it{NL}','fontsize',4);
m_text(-63.0546-0.3,49.4768,'\it{AI}','fontsize',6);
m_text(-66.2018,48.6351,'\it{QC}','fontsize',4);
m_text(-59.2638+0.4,47,{'CS'},'fontsize',4,...
    'FontAngle', 'italic','HorizontalAlignment','center','color','k');
m_text(-59.1171,50.5,'\it{BI}','fontsize',4);
m_annotation('textarrow',[-68.95 -69.6811],[49.25 48.1688+0.1],...
    'string','\it{TS}','fontsize',6,'Fontweight','normal',...
    'horizontalalignment','center','headlength',5);
m_annotation('textarrow',[-68.95+1.3 -69.6811+1.5],[49.25+0.5 48.1688+0.85],...
    'string','\it{LSLE}','fontsize',6,'Fontweight','normal',...
    'horizontalalignment','center','headlength',5);
m_scatter(-69.681130,48.168871,10,'filled','k');
m_text(-69.3,47,'a)','fontsize',8);
s(1)=m_scatter(NaN,NaN,15,'^','filled',...
    'markeredgecolor','none','markerfacecolor','k');
s(2)=m_scatter(NaN,NaN,15,'o','filled',...
    'markeredgecolor','none','markerfacecolor','k');

set(findall(gcf,'-property','fontname'),'fontname','Latin Modern Roman');
set(findall(gcf,'-property','fontsize'),'fontsize',8);

ax2=axes('position',[0.65 0.25 0.15 0.7]); hold on;
fill([GoSLmnTemp-2*GoSLstdTemp;flipud(GoSLmnTemp+2*GoSLstdTemp)],...
    [1:400 400:-1:1],rgb_x('deep pink'),'facealpha',0.1,'edgecolor','none');
plot(CTD(17).t090C,CTD(17).depSM*-1,'color',rgb_x('deep pink'),'linewidth',1);
ax2.YAxis.Visible='off';ax2.Color='none'; ax2.XAxis.Color=rgb_x('deep pink'); 
axis ij;axis tight; ylim([0 400]);
xticks([0;5;10]); xticklabels({'0';'5';'10 °C'});
xtickangle(0);
text(1,375,'c)');

ax3=axes('position',[0.65+0.05 0.25 0.15 0.7]); hold on;
fill([GoSLmnSR-2*GoSLstdSR;flipud(GoSLmnSR+2*GoSLstdSR)],...
    [1:400 400:-1:1],rgb_x('cornflower'),'facealpha',0.1,'edgecolor','none');
plot(CTD(17).sr,CTD(17).depSM*-1,'color',rgb_x('cornflower'),'linewidth',1);
axis ij;axis tight; ylim([0 400]);xticks(30:2.5:35);
xt=ax3.XAxis.TickValues;xl=xlim;
ax3.YAxis.Visible='off'; 
ax3.Color='none'; 
ax3.XAxis.Visible='off'; 

iax3=axes('position',[0.65+0.05 0.25-0.075 0.15 0.7]);
iax3.YAxis.Visible='off'; 
iax3.Color='none'; 
iax3.XAxis.Color=rgb_x('cornflower');
xlim(xl);xticks(xt);
xticklabels({'30';'32.5';'35 g/kg'});
xtickangle(0);

ax4=axes('position',[0.65+0.1 0.25 0.15 0.7]); hold on;
fill([GoSLmnDO-2*GoSLstdDO;flipud(GoSLmnDO+2*GoSLstdDO)],...
    [1:400 400:-1:1],rgb_x('charcoal'),'facealpha',0.1,'edgecolor','none');
plot(CTD(17).sbox0Mm,CTD(17).depSM*-1,'color',rgb_x('charcoal'),'linewidth',1);
ax4.YAxis.Visible='off';
ax4.Color='none'; 
axis ij;axis tight; ylim([0 400]);xticks(150:100:350);
xt=ax4.XAxis.TickValues;xl=xlim;
ax4.YAxis.Visible='off'; 
ax4.Color='none'; 
ax4.XAxis.Visible='off'; 

iax4=axes('position',[0.65+0.1 0.25-0.15 0.15 0.7]);
iax4.YAxis.Visible='off'; 
iax4.Color='none'; 
iax4.XAxis.Color=rgb_x('charcoal');
xlim(xl);xticks(xt);
xticklabels({'150';'250';'      350 \mumol/kg'});
xtickangle(0);

iax=axes('Position',[0.65 0.25 0.30 0.7]);
linkaxes([iax ax2 ax3 ax4],'y'); 
axis ij;
iax.XAxis.Visible='off';
uistack(iax,'bottom');
iax.Color='w'; iax.YGrid='on';
iax.GridAlpha=0.05;
box on
yticks([0:100:400]);
box off; xl=xlim;
line([0 xl(2)],[275 275],'color','k','linestyle','--','linewidth',0.5);
ylabel('\it{z} \rm{(m)}');

h1=histcounts(vertcat(CTD1.zStar),0:10:500);
h2=histcounts(vertcat(CTD2.zStar),0:10:500);
ax5=axes('position',[0.65+0.25 0.25 0.05 0.7]);
b(1)=barh(0+5:10:500-5,h1,1,'facecolor','b','facealpha',0.3,'edgecolor','none');
hold on
b(2)=barh(0+5:10:500-5,h2,1,'facecolor','r','facealpha',0.3,'edgecolor','none');
axis ij tight; ylim([0 400]); xticks(0:25:50);

ax5.Color='none'; 
% ax5.YAxis.Visible='on';
yticklabels([]);
ax5.Box='off';
xlabel('\it{n}');
ylabel('\it{z^*} \rm{(m)}');
set(findall(gcf,'-property','fontsize'),'fontsize',8);
scatter(30,50,20,'s','markerfacecolor','b','markerfacealpha',0.3,...
    'markeredgecolor','none');
text(39,49,'Leg 1','fontsize',6);
scatter(30,65,20,'s','markerfacecolor','r','markerfacealpha',0.3,...
    'markeredgecolor','none');
text(39,64,'Leg 2','fontsize',6);
text(8,375,'d)','fontsize',8);

axG=axes('position',[.06 .725 .25 .25]); hold on;
m_proj('lambert','long',[-75 -45],'lat',[40 55]);
[CS,CH]=m_etopo2('contour',[-5000 -2500 -100:100:0],'color',[0.9 0.9 0.9]);
m_gshhs_i('patch',col,'edgecolor',[0.3 0.3 0.3],'linewidth',0.3);
m_grid('linest','-','xtick',[],'ytick',[],'linewidth',...
    0.2,'fontsize',5);
m_scatter(dLon,dLat,50,'p','filled',...
    'markerfacecolor',rgb_x('light violet'),'markeredgecolor','k'); 
m_rectangle(lon_lim(1),lat_lim(1),diff(lon_lim),diff(lat_lim),...
    0,'edgecolor',rgb_x('black'),'linewidth',0.5);
m_text(-68,42,'b)','fontsize',8);
m_text(-52,42.5,{'\it{Atlantic}';'\it{Ocean}'},'fontsize',5,'horizontalalignment','center');

axes(ax1)
legend(s,'Interior    (\it{h} \rm{>279 m)}','Boundary (\it{h} \rm{<279 m)}',...
    'location','southwest','fontsize',6);

axes(ax4);
scatter(300,300,20,'s','markerfacecolor',rgb_x('deep pink'),...
    'markerfacealpha',1,'markeredgecolor','none');
text(310,300,'Temp.','fontsize',6);

scatter(300,315,20,'s','markerfacecolor',rgb_x('cornflower'),...
    'markerfacealpha',1,'markeredgecolor','none');
text(310,315,'Sal.','fontsize',6);

scatter(300,330,20,'s','markerfacecolor',rgb_x('charcoal'),...
    'markerfacealpha',1,'markeredgecolor','none');
text(310,330,'D.O.','fontsize',6);
axes(axG)

set(findall(gcf,'-property','fontname'),'fontname','Latin Modern Roman');
%% 
export_fig /Users/samst/Dropbox/UBC/TReX_paper/figures/R1/map_R1.pdf -dpdf -nofontswap