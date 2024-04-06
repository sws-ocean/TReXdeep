%% Along channel distribution
addpath(genpath('/Users/samst/Dropbox/UBC/TReX_deepII/'));
addpath(genpath('/Users/samst/Dropbox/UBC/TReX_paper/'));
addpath(genpath('/Users/samst/Dropbox/UBC/GSW/'));
addpath(genpath('/Users/samst/Dropbox/UBC/m_map/'));
addpath(genpath('/Users/samst/Dropbox/UBC/Misc/'));

clear 

bathy=load('GoSL_Zhires.mat');
bathy.Z(bathy.Z==-9999)=NaN;

dec=5;
bathy.lon=bathy.lon(1:dec:end,1:dec:end);
bathy.lat=bathy.lat(1:dec:end,1:dec:end);
bathy.Z=bathy.Z(1:dec:end,1:dec:end);

load data20230707.mat

%% Find all stations in Laurentian channel
lat_lim=[46 51];
lon_lim=[-70 -57.5];
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

LC.lon=[-66.4852,-65.8573,-65.1713,-64.4853,-63.9970,-63.4506,-64.2761,...
    -65.4969,-66.1712,-66.9735,-63.4622,-63.9738,-63.8110,-62.8460,...
    -62.5669,-63.4506,-62.2530,-62.5786,-61.9507,-61.4275,-60.6601,...
    -61.4159,-62.5786,-62.8576,-60.6601,-61.4159,-60.7996,-60.1602,...
    -59.3463,-59.3928,-59.6253,-57.5674,-58.7998,-58.0441,-59.8579,...
    -59.7416-0.2,-58.7766,-60.1020,-60.4160,-60.9392,-60.6601,...
    -59.6253,-59.4858,-59.4044,-59.1719,-60.1485,-59.7893];

LC.lat=[49.3725,49.2800,49.3032,49.1568,48.8794,49.3802,49.7808,49.8964,...
    50.0813,49.7038,49.3802,48.8794,48.5636,48.4711,49.1183,49.3725,...
    49.1106+0.1,48.2246+0.1,48.2246,48.0859,48.7870,...
    48.9488,49.1106,48.4634,48.7870,48.0859,47.8163,47.2461,47.5851,...
    48.1629,48.3401,46.5913,46.0289,47.2153,46.6529,49.6653-0.1,...
    49.5343-0.25,49.2646,49.1029,48.8717,48.7870,...
    48.3401,48.2631,48.5867,49.0412,49.2569,48.8748];

% Lazy point removal
pl=[28;29;32;33;34;35;36;37;38;39;44;45;46];
LC.lon(pl)=[];LC.lat(pl)=[];

LC.station_bound=alphaShape(LC.lon',LC.lat',1,'HoleThreshold',15);
[l,ll]=m_ll2xy(LC.lon,LC.lat);
tmpshp=alphaShape(l',ll',0.01);
LC.AS=plot(tmpshp,'facealpha',0.1,'edgecolor','none','facecolor',rgb_x('forest green'));

idx1=find(inShape(LC.station_bound,stlon1,stlat1));
idx2=find(inShape(LC.station_bound,stlon2,stlat2));

% Pull out channel concs and distance
ICTR1=sum(vertcat(CTD1(idx1).I_j).*1026,2);
ICTR2=sum(vertcat(CTD2(idx2).I_j).*1026,2);    
ICx1=vertcat(CTD1(idx1).x_dist); 
ICx2=vertcat(CTD2(idx2).x_dist); 

stdist1=unique(ICx1);
mnTR1=NaN(length(stdist1),1);
for i=1:length(stdist1)
    mnTR1(i)=ICTR1(ICx1==stdist1(i));
end

stdist2=unique(ICx2);
mnTR2=NaN(length(stdist2),1);
% Average together core samples
for i=1:length(stdist2)
    mnTR2(i)=ICTR2(ICx2==stdist2(i));
end

% Make in-channel bin average
inter1=50;
channelX1=0:inter1/2:700;
centredX1=inter1/2:inter1/2:max(channelX1)+(inter1/2);
count=0;
ICmn1=NaN(length(channelX1),1);
ICstd1=NaN(length(channelX1),1);ICstd2=NaN(length(channelX1),1);
for i=channelX1
    count=count+1;
    ICmn1(count)=nanmean(mnTR1(stdist1>i & stdist1<i+inter1));
    ICstd1(count)=std(mnTR1(stdist1>i & stdist1<i+inter1),'omitnan');
end

inter2=75;
channelX2=0:inter2/2:700;
centredX2=inter2/2:inter2/2:max(channelX2)+(inter2/2);
count=0;
ICmn2=NaN(length(channelX2),1);
ICstd2=NaN(length(channelX2),1);
for i=channelX2
    count=count+1;
    ICmn2(count)=nanmean(mnTR2(stdist2>i & stdist2<i+inter2));
    ICstd2(count)=std(mnTR2(stdist2>i & stdist2<i+inter2),'omitnan');
end

figure;
axes;hold on;
scatter(stdist1,mnTR1,20,'filled','markerfacecolor','b','markeredgecolor',...
    'none','markerfacealpha',0.3);
plot(centredX1,ICmn1,'b','linewidth',1.5);
scatter(centredX1,ICmn1,30,'filled','markerfacecolor','b',...
    'markeredgecolor','k','markerfacealpha',0.3);
% errorbar(centredX1,ICmn1,ICstd1,'markersize',5,'marker','o',...
%     'markerfacecolor','b','markeredgecolor','k',...
%     'color','b');

scatter(stdist2,mnTR2,20,'filled','markerfacecolor','r','markeredgecolor',...
    'none','markerfacealpha',0.3);
plot(centredX2,ICmn2,'r','linewidth',1.5);
scatter(centredX2,ICmn2,30,'filled','markerfacecolor','r',...
    'markeredgecolor','k','markerfacealpha',0.3);
% errorbar(centredX2,ICmn2,ICstd2,'markersize',5,'marker','o',...
%     'markerfacecolor','r','markeredgecolor','k',...
%     'color','r');


%% Do the same for 3d model data
mTracer=load('conc_20220604_K10');
mTracer.x_dist=gsw_distance([ones(size(mTracer.lon(:)))*dLon mTracer.lon(:)]...
    ,[ones(size(mTracer.lat(:)))*dLat mTracer.lat(:)])/1000;
idx3=find(inShape(LC.station_bound,mTracer.lon,mTracer.lat));

mChannel=mTracer.C(idx3);
mDist=mTracer.x_dist(idx3);

mDist(isnan(mChannel))=[];
mChannel(isnan(mChannel))=[];

% Make in-channel bin average
interM=50;
channelXM=0:interM/2:700;
centredXM=interM/2:interM/2:max(channelXM)+(interM/2);
count=0;
ICmnM=NaN(length(channelXM),1);
ICstdM=NaN(length(channelXM),1);
for i=channelXM
    count=count+1;
    ICmnM(count)=nanmean(mChannel(mDist>i & mDist<i+interM));
    ICstdM(count)=std(mChannel(mDist>i & mDist<i+interM),'omitnan');
end

% %% How much tracer is in the LC?
% idx4=find(inShape(LC.station_bound,bathy.lon,bathy.lat));
% allZ=bathy.Z(idx4)*-1;
% allA=cdtarea(bathy.lat,bathy.lon);allA=allA(idx4);
% Az=zeros(ceil(abs(max(bathy.Z(:)))),1);
% Az_z=0:max(allZ);
% for i=0:max(allZ)
%     Az(i+1)=sum(allA(allZ>i));
% end
% 
% tmp1=horzcat(CTD1.c_ij);
% tmp2=horzcat(CTD2.c_ij);
% 
% Azi=interp1(Az_z'-275,Az,zHat_grid);
% MLC1=simps(zHat_grid,Azi.*mean(tmp1(:,idx1),1:2),1)*1026.25/1e15*mm_sf5cf3
% MLC2=simps(zHat_grid,Azi.*mean(tmp2(:,idx2),1:2),1)*1026.25/1e15*mm_sf5cf3


%% can I model? yes!
% THese results currently (A=126 m2/s, U= 400 m/day) are similar to our
% estimates of A in SoG IW and to landward along-channel flow estimates
% (see Gilbert 2005, Bugden 1988). 

M0=[71 60]; % g
mm_sf5cf3 = 196.06;
X=0:1e3:725e3; % m
A=75:1:250; %m2/s
U=100/86400:30/86400:750/86400;% m/s
t=[(CTD1(1).mtime(1)-datenum(2021,10,24))*86400;...
    (CTD2(1).mtime(1)-datenum(2021,10,24))*86400]; % s

ychannel=75e3; zspread=1; % m
Achannel=ychannel;%*zspread; % m2, area of the channel*depth range of tracer
c1=NaN(length(X),length(U)*length(A));cost1=NaN(1,length(U)*length(A));
c2=NaN(length(X),length(U)*length(A));cost2=NaN(1,length(U)*length(A));
cm=NaN(length(X),length(U)*length(A));costM=NaN(1,length(U)*length(A));
idxs=NaN(2,length(U)*length(A));
count=0;

[tmpDist1,idxD1]=sort(stdist1);
tmpTR1=mnTR1;tmpTR1=tmpTR1(idxD1);
[tmpDist2,idxD2]=sort(stdist2);
tmpTR2=mnTR2;tmpTR2=tmpTR2(idxD2);
[tmpDistM,idxM1]=sort(mDist);
tmpTRM=mChannel;tmpTRM=tmpTRM(idxM1);

cost1N=NaN(1,length(U)*length(A));
cost2N=NaN(1,length(U)*length(A));

ICmn1(isnan(ICmn1))=0;
ICmn2(isnan(ICmn2))=0;

for i=1:length(A)
    for ii=1:length(U)
        count=count+1;
        c1(:,count)=M0(1)./sqrt(4*pi*A(i)*t(1)).*exp( -(X-U(ii)*t(1)).^2./(4*A(i)*t(1)))...
            ./1e3.*1e15./mm_sf5cf3/Achannel.*1026; % fmol/m2
        cm(:,count)=1./sqrt(4*pi*A(i)*t(1)).*exp( -(X-U(ii)*t(1)).^2./(4*A(i)*t(1)))...
            *1e6/Achannel;
        c2(:,count)=M0(2)./sqrt(4*pi*A(i)*t(2)).*exp( -(X-U(ii)*t(2)).^2./(4*A(i)*t(2)))...
            ./1e3.*1e15./mm_sf5cf3/Achannel.*1026; % fmol/m2
        cost1(count)=trapz(tmpDist1,(interp1(X,c1(:,count),tmpDist1*1e3)-...
            tmpTR1).^2);
        cost2(count)=trapz(tmpDist2,(interp1(X,c2(:,count),tmpDist2*1e3)-...
            tmpTR2).^2);
        costM(count)=trapz(tmpDistM,(interp1(X,cm(:,count),tmpDistM*1e3)-...
            tmpTRM).^2);
        idxs(:,count)=[i;ii];
        
        c1N=(c1(:,count)-min(c1(:,count)))/(max(c1(:,count))-...
            min(c1(:,count))).*max(ICmn1);
        cost1N(count)=trapz(centredX1*1000,(interp1(X,c1(:,count),centredX1*1000)-...
            ICmn1').^2);

        c2N=(c2(:,count)-min(c2(:,count)))/(max(c2(:,count))-...
            min(c2(:,count))).*max(ICmn2);
        cost2N(count)=trapz(centredX2*1000,(interp1(X,c2(:,count),centredX2*1000)-...
            ICmn2').^2);

    end
end
    
[~,idxC1]=min(cost1); [~,idxC2]=min(cost2); [~,idxCM]=min(costM);
[~,idxC1N]=min(cost1N); [~,idxC2N]=min(cost2N);

c2from1=M0(1)./sqrt(4*pi*A(idxs(1,idxC1))*t(2)).*...
    exp( -(X-U(idxs(2,idxC1))*t(2)).^2./(4*A(idxs(1,idxC1))*t(2)))...
    ./1e3.*1e15./mm_sf5cf3/Achannel; % fmol/l

fprintf('\nLeg 1 1D model: k=%3.0f m^2 s^{-1}, U=%2.2f cm s^{-1}\n',...
    A(idxs(1,idxC1)),U(idxs(2,idxC1))*100);
fprintf('\nLeg 1 1D model (CIOPS-E): k=%3.0f m^2 s^{-1}, U=%2.2f cm s^{-1}\n',...
    A(idxs(1,idxCM)),U(idxs(2,idxCM))*100);
fprintf('\nLeg 2 1D model: k=%3.0f m^2 s^{-1}, U=%2.2f cm s^{-1}\n',...
    A(idxs(1,idxC2)),U(idxs(2,idxC2))*100);

fprintf('\nLeg 1 NORMALIZED 1D model: k=%3.0f m^2 s^{-1}, U=%2.2f cm s^{-1}\n',...
    A(idxs(1,idxC1N)),U(idxs(2,idxC1N))*100);

fprintf('\nLeg 2 NORMALIZED 1D model: k=%3.0f m^2 s^{-1}, U=%2.2f cm s^{-1}\n',...
    A(idxs(1,idxC2N)),U(idxs(2,idxC2N))*100);


cf1=c1(:,idxC1);cf2=c2(:,idxC1);
save('/Users/samst/Library/CloudStorage/Dropbox/UBC/TReX_paper/code/letter/R1/channel1D_R1.mat',...
    'centredX1','ICmn1','centredX2','ICmn2','X','cf1','cf2');

%% Bootstrap some errors on the model 
I1=bootrnd(size(tmpTR1,1),1000);
I1M=bootrnd(size(tmpTRM,1),1000);
I2=bootrnd(size(tmpTR2,1),1000);

% % I think this method is a Jackknife, removes one sample per series
% I1=randi(size(tmpTR1,1),size(tmpTR1,1),1000);
% I1(1,:)=[];
% I2=randi(size(tmpTR2,1),size(tmpTR2,1),1000);
% I2(1,:)=[];

Aerr1=NaN(size(I1,2),1); Uerr1=NaN(size(I1,2),1);
Aerr1M=NaN(size(I1M,2),1); Uerr1M=NaN(size(I1M,2),1);
Aerr2=NaN(size(I2,2),1); Uerr2=NaN(size(I2,2),1);

hw=waitbar(0,'Calculating...');
tic
for j=1:size(I1,2)
    idxsE=NaN(2,length(U)*length(A));
    count=0;
    cost1E=NaN(1,length(U)*length(A));
    cost1EM=NaN(1,length(U)*length(A));
    cost2E=NaN(1,length(U)*length(A));
    c1E=NaN(length(X),length(U)*length(A));
    c1EM=NaN(length(X),length(U)*length(A));
    c2E=NaN(length(X),length(U)*length(A));
    for i=1:length(A)
        for ii=1:length(U)
            count=count+1;
            c1E(:,count)=M0(1)./sqrt(4*pi*A(i)*t(1)).*exp( -(X-U(ii)*t(1)).^2./(4*A(i)*t(1)))...
                ./1e3.*1e15./mm_sf5cf3/Achannel.*1026; % fmol/m2
            c1EM(:,count)=1.2./sqrt(4*pi*A(i)*t(1)).*exp( -(X-U(ii)*t(1)).^2./(4*A(i)*t(1)))...
                *1e6/Achannel;
            c2E(:,count)=M0(2)./sqrt(4*pi*A(i)*t(2)).*exp( -(X-U(ii)*t(2)).^2./(4*A(i)*t(2)))...
                ./1e3.*1e15./mm_sf5cf3/Achannel.*1026; % fmol/m2
            
            D1=tmpDist1(I1(:,j));
            T1=tmpTR1(I1(:,j));
            [D1,idxD1]=sort(D1);
            T1=T1(idxD1);
            cost1E(count)=trapz(D1,(interp1(X,c1E(:,count),D1*1e3)-...
                T1).^2);
            
            D1=tmpDistM(I1M(:,j));
            T1=tmpTRM(I1M(:,j));
            [D1,idxD1]=sort(D1);
            T1=T1(idxD1);
            cost1EM(count)=trapz(D1,(interp1(X,c1EM(:,count),D1*1e3)-...
                T1).^2);
            
            D2=tmpDist2(I2(:,j));
            T2=tmpTR2(I2(:,j));
            [D2,idxD2]=sort(D2);
            T2=T2(idxD2);
            cost2E(count)=trapz(D2,(interp1(X,c2E(:,count),D2*1e3)-...
                T2).^2);
            
            idxsE(:,count)=[i;ii];
        end
    end 
    [~,idxC1E]=min(cost1E); [~,idxC2E]=min(cost2E);[~,idxC1EM]=min(cost1EM);
    
    Aerr1(j,1)=A(idxsE(1,idxC1E));
    Uerr1(j,1)=U(idxsE(2,idxC1E));
    Aerr1M(j,1)=A(idxsE(1,idxC1EM));
    Uerr1M(j,1)=U(idxsE(2,idxC1EM));
    Aerr2(j,1)=A(idxsE(1,idxC2E));
    Uerr2(j,1)=U(idxsE(2,idxC2E));
       
    waitbar(j/size(I1,2),hw);
end
toc
delete(hw)           
% 
% 
% save('/Users/samst/Library/CloudStorage/Dropbox/UBC/TReX_paper/data/1Dboot.mat',...
%     'Aerr1', 'Uerr1', 'Aerr1M', 'Uerr1M', 'Aerr2', 'Uerr2');
% % 
% c1Elim(:,1)=M0(1)./sqrt(4*pi*(A(idxs(1,idxC1))-Aerr1)*t(1)).*...
%     exp( -(X-(U(idxs(2,idxC1))-Uerr1)...
%     *t(1)).^2./(4*(A(idxs(1,idxC1))-Aerr1)*t(1)))...
%     ./1e3.*1e15./mm_sf5cf3/Achannel; % fmol/m2
% 
% c1Elim(:,2)=M0(1)./sqrt(4*pi*(A(idxs(1,idxC1))+Aerr1)*t(1)).*...
%     exp( -(X-(U(idxs(2,idxC1))+Uerr1)...
%     *t(1)).^2./(4*(A(idxs(1,idxC1))+Aerr1)*t(1)))...
%     ./1e3.*1e15./mm_sf5cf3/Achannel; % fmol/m2
% 
% c1Elim(:,3)=M0(1)./sqrt(4*pi*(A(idxs(1,idxC1))-Aerr1)*t(1)).*...
%     exp( -(X-(U(idxs(2,idxC1))+Uerr1)...
%     *t(1)).^2./(4*(A(idxs(1,idxC1))-Aerr1)*t(1)))...
%     ./1e3.*1e15./mm_sf5cf3/Achannel; % fmol/m2
% 
% c1Elim(:,4)=M0(1)./sqrt(4*pi*(A(idxs(1,idxC1))+Aerr1)*t(1)).*...
%     exp( -(X-(U(idxs(2,idxC1))-Uerr1)...
%     *t(1)).^2./(4*(A(idxs(1,idxC1))+Aerr1)*t(1)))...
%     ./1e3.*1e15./mm_sf5cf3/Achannel; % fmol/m2
% 
% c1ElimF(:,1)=min(c1Elim,[],2);
% c1ElimF(:,2)=max(c1Elim,[],2);
% 
% l(2)=plot(X/1e3,c1(:,idxC1),'b','linewidth',2);
% hold on; plot(X/1000,c1ElimF,'k');

load 1Dboot.mat

% From 1000 subsets
fprintf('\nLeg 1 (2 STD): k=%3.0f m^2 s^{-1}, U=%4.3f cm s^{-1}\n',...
    std(Aerr1)*2,std(Uerr1)*2*100);
fprintf('\nLeg 1 CIOPSE (2 STD): k=%3.0f m^2 s^{-1}, U=%4.3f cm s^{-1}\n',...
    std(Aerr1M)*2,2*std(Uerr1M)*100);
fprintf('\nLeg 2 (2 STD): k=%3.0f m^2 s^{-1}, U=%4.3f cm s^{-1}\n',...
    std(Aerr2)*2,std(Uerr2)*2*100);

%% Pretty plot
f1=figure('units','centimeters','outerposition',...
     [0.01 0.01 20 15],'color','w');
lat_lim=[47 51];
lon_lim=[-70 -57.5];
col=[255 214 140]/255; % YELLOW!

lnes=linspecer(7); lnes([1:2],:)=[];

%%%%%%% Axis and stations %%%%%%%
ax1=axes('position',[.05 .55 0.5 0.4]); hold on;
m_proj('oblique','lon',fliplr(lon_lim),'lat',lat_lim,'dir','horz','aspect',0.6)
[c,h]=m_contour(bathy.lon,bathy.lat,bathy.Z,0:-100:-500,'color',...
    [0.95 0.95 0.95],'linewidth',0.5);
clabel(c,h,'LabelSpacing',500,'color',[0.95 0.95 0.95],'fontsize',6);
m_gshhs_f('patch',col,'edgecolor',[0.3 0.3 0.3]); % CHANGE
m_grid('linestyle','none','linewidth',0.1,'tickdir','out',...
    'xaxisloc','top','yaxisloc','left','fontsize',6);
LC.AS=plot(tmpshp,'facealpha',0.05,'edgecolor','none','facecolor',...
    rgb_x('lilac'));
m_scatter(stlon1(idx1),stlat1(idx1),25,'filled',...
    'markerfacecolor','b','markeredgecolor','w','markerfacealpha',0.15);
m_scatter(stlon2(idx2),stlat2(idx2),25,'filled',...
    'markerfacecolor','r','markeredgecolor','w','markerfacealpha',0.15);

text(1.9293,0.0406,'Leg 1','fontsize',6,'color',rgb_x('cornflower'));
text(1.9293,0.035,50.8,'Leg 2','fontsize',6,'color',rgb_x('deep pink'));


% Build axis
x=linspace( -59.9288,-65.9100,1000);
y=linspace(48.2310,49.8008,1000);
d=gsw_distance([ones(1,length(x))*-59.4858;x]',...
    [ones(1,length(y))*48.7254;y]')/1000;
[~,idx]=min(abs(d-500));
x(idx:end)=[]; y(idx:end)=[];d(idx:end)=[];
x2=linspace(x(end),-69.4518);
y2=linspace(y(end),48.1003);
d2=gsw_distance([ones(1,length(x2))*x(end);x2]',...
    [ones(1,length(y2))*y(end);y2]')/1000;
x=[x x2];y=[y y2]; d=[d;d(end)+d2];

m_plot(x,y,'k','linewidth',1.5);

% m_scatter(-69.670562,48.168191,50,'filled')

count=0;
for i=50:50:max(d)
    [~,idx]=min(abs(d-i));
    
    if rem(i,200)==0
        count=count+1;
        m_scatter(x(idx),y(idx),50,'filled','markerfacecolor',lnes(count,:),...
        'markeredgecolor','k','marker','s');
    else
        m_scatter(x(idx),y(idx),10,'filled','markerfacecolor','k',...
            'markeredgecolor','k');
    end
end
%     if rem(i,200)==0 && i<500
%         count=count+1;
%         m_text(x(idx),y(idx),num2str(i),'horizontalalignment','center',...
%             'rotation',10,'color',lnes(count,:),'backgroundcolor','w');
%     elseif rem(i,200)==0 && i>500
%         count=count+1;
%         m_text(x(idx),y(idx),num2str(i),'horizontalalignment','center',...
%             'rotation',0,'color',lnes(count,:),'backgroundcolor','w');
%     end

s=m_scatter(dLon,dLat,200,'p','filled',...
    'markerfacecolor',[0.929,0.694,0.125],'markeredgecolor','k');
m_northarrow(-66.4330,51.5,1,'type',3,'linewi',.5);
m_text(-67.8003,53.0036,'a)','fontsize',8);

%%%%%%% spreading %%%%%%%
l=[];
ax2=axes('position',[.1 .2 0.4 0.3]); hold on;
s(1)=scatter(stdist1,mnTR1,20,'filled','markerfacecolor','b','markeredgecolor',...
    'none','markerfacealpha',0.175);
plot(centredX1,ICmn1,'b--','linewidth',1);
scatter(centredX1,ICmn1,30,'marker','s','markerfacecolor','b',...
    'markeredgecolor','none','markerfacealpha',0.5);
s(2)=scatter(stdist2,mnTR2,20,'filled','markerfacecolor','r','markeredgecolor',...
    'none','markerfacealpha',0.175);
plot([stdist2(1) centredX2],[mnTR2(1);ICmn2],'r--','linewidth',1);
scatter([stdist2(1) centredX2],[mnTR2(1);ICmn2],30,'marker','s','markerfacecolor','r',...
    'markeredgecolor','none','markerfacealpha',0.6);
lp(1)=scatter(NaN,NaN,20,'markerfacecolor','k',...
    'markeredgecolor','none','markerfacealpha',0.175);
lp(2)=plot(NaN,NaN,'k--','marker','s','markerfacecolor','k',...
    'markeredgecolor','none');
% plot(X/1e3,c1(:,idxC1),'w','linewidth',2);
l(2)=plot(X/1e3,c1(:,idxC1),'b-','linewidth',1.5);
lp(3)=plot(NaN,NaN,'k-','linewidth',1.5);
% plot(X/1e3,c2(:,idxC2),'w','linewidth',2);
l(4)=plot(X/1e3,c2(:,idxC2),'r-','linewidth',1.5);
set(gca,'ycolor','k');
xlabel('\it{x} \rm (km)');
ylabel('\it{I_j} \rm(fmol m^{-2})');
axis tight; grid on;%set(gca,'xdir','reverse');
xlim([0 525]);
text(15,37000,'b)','fontsize',8);

set(findall(gcf,'-property','fontsize'),'fontsize',8);
set(findall(gcf,'-property','fontname'),'fontname','Latin Modern Roman');
 


% legend(l,'Overlapping window mean','location','best');

% ax1=axes('position',[.2 .5 1 1]); hold on;


%%%%%%%% SO how long is it going to take to get up the estuary %%%%%%%

y=[linspace(-59.9353,-66.0647,200) linspace(-66.0647-0.01,-68.8538,200)];
yy=[linspace(48.2394,49.5798,200) linspace(49.5798-0.01,48.4947,200)];
transX=cumsum(gsw_distance(y,yy));

l=[-61.7190,-64.5080,-67.199,-68.8538];
ll=[48.6011,49.3883,49.1755,48.4947];

X1d=0:1e3:800e3; % m
A_m=mean([A(idxs(1,idxC1)) A(idxs(1,idxC2))]); %m2/s
U_m=mean([U(idxs(2,idxC1)) U(idxs(2,idxC1))]);% m/s
t=86400:86400:15*365*86400; % s

ychannel=75e3; zspread=1; % m
Achannel=ychannel;%*zspread; % m2, area of the channel*depth range of tracer

cm1=NaN(length(t),length(X1d));
for i=1:length(t)
    cm1(i,:)=M0(1)./sqrt(4*pi*A_m*t(i)).*exp(-(X1d-U_m*t(i)).^2./(4*A_m*t(i)))...
        ./1e3.*1e15./mm_sf5cf3/Achannel*1026; % fmol/m2
end

% dd=interp1(y(2:end), transX,l);dd(end)=max(transX);

ax3=axes('position',[.545 .2 0.15 0.7]); 
[xx,yy]=meshgrid(X1d,t);
contourf(xx./1e3,yy./86400./365,cm1,[0:0.05:0.5 1e3:0.5e3:20e3],'linestyle','none');
% caxis([1 20]);
tmp=seminfhaxby;tmp2=linspecer;%.*repmat(1-.2*cos([0:127]'/64*2*pi).^4,1,3);
colormap(gca,flipud([flipud(tmp2(1:end,:));tmp(1,:)]));
cc=colorbar;
xlabel('\it{x} \rm (km)');
ylabel('\it{t} \rm(years)');
% title(cc,'SF5 (fmol m^{-2})');
text(50,9.7,'c)','fontsize',8);
ylim([0 10]);
xticks(0:250:750);

ax4=axes('position',[.74 .2 0.18 0.7]);
hold on

dd=[200:200:800];

for i=1:length(dd)
    %     [~,idx]=min(abs(X-dd(i)));
    [~,idx]=min(abs(X/1000-dd(i)));
    % /1000 unit conversion added
    plot(t/86400/365,cm1(:,idx)/10000,'color',lnes(i,:),'linewidth',2);
end
grid on; axis tight;
% le=legend('1','2','3','4','location','best');
xlabel('\it{t} \rm(years)');
ylabel('\it{I''_j} \rm(fmol m^{-2})');
% title(le,'Station')
text(2.2,1.7,'d)','fontsize',8);
xlim([0 10]);
xticks(0:2:10);
text(-0.5,1.8,'×10^4');

set(findall(gcf,'-property','fontsize'),'fontsize',8);
axes(ax2)
legend(lp,'Observations','Window mean','1D model','location','northeast',...
    'fontsize',6);
set(findall(gcf,'-property','fontname'),'fontname','Latin Modern Roman');

%%
export_fig /Users/samst/Dropbox/UBC/TReX_paper/figures/R1/channel_thesisR1.pdf -dpdf -nofontswap
% 
%%
% 
% axes;
% hold on;
% lat_lim=[46 52];
% lon_lim=[-70 -56];
% m_proj('equidistant','lon',lon_lim,'lat',lat_lim);
% [CS,CH]=m_etopo2('contourf',0:-100:-600,'linestyle','none');
% caxis([-500 0]);
% colormap(cm);
% m_gshhs_h('patch',col,'edgecolor','k');
% m_grid('linestyle','none','linewidth',0.1,'tickdir','out',...
%     'xaxisloc','top','yaxisloc','right','fontsize',6);
% s=m_scatter(dLon,dLat,200,'p','filled',...
%     'markerfacecolor',[0.929,0.694,0.125],'markeredgecolor','k');
% 
% m_plot(y,yy,'w');
% 
% for i=1:4
%     m_scatter(l(i),ll(i),100,'filled','markerfacecolor',lnes(i,:),...
%         'markeredgecolor','w','markerfacealpha',1);
%     m_text(l(i)+0.2,ll(i),num2str(i),'color','w',...
%         'fontweight','bold');
% end
% 
% % m_text(-58.0872-0.5,48.4215-0.4,'\itNewfoundland','fontsize',4);
% % % m_text(-63.0546-0.3,49.4768,'\it{Anticosti Island}','fontsize',6);
% % m_text(-66.4483-1.25,47.5129-.5,'\it{New Brunswick}','fontsize',4);
% % % m_text(-63.5096,46.3486,'\it{P.E.I.}','fontsize',6);
% % m_text(-66.2018,48.6351,'\it{Québec}','fontsize',4);
% % m_text(-60,49.458724,{'Gulf of';'St. Lawrence'},'fontsize',4,...
% %     'FontAngle', 'italic','HorizontalAlignment','center','color','w');
% % m_text(-59.2638,47.1793,{'Cabot';'Strait'},'fontsize',4,...
% %     'FontAngle', 'italic','HorizontalAlignment','center','color','w');
% 
% 
% set(findall(gcf,'-property','fontsize'),'fontsize',12);
% set(findall(gcf,'-property','fontname'),'fontname','Latin Modern Roman');
% %%
% % export_fig /Users/samst/Dropbox/UBC/TReX_paper/figures/channelTrans.pdf -pdf -nofontswap
% 
% %% 
% 
% figure('color','w');
% 
% 
% 
% subplot(1,2,2); hold on
% for i=1:length(dd)
%     [~,idx]=min(abs(X-dd(i)));
%     plot(t/86400/365,cm1(:,idx),'color',lnes(i,:),'linewidth',2);
% end
% grid on; axis tight;
% le=legend('1','2','3','4','location','best');
% xlabel('\it{t} \rm(years)');
% ylabel('SF5 (fmol m^{-2})');
% title(le,'Station')
% 
% set(findall(gcf,'-property','fontsize'),'fontsize',12);
% set(findall(gcf,'-property','fontname'),'fontname','Latin Modern Roman');
% %%
% % export_fig /Users/samst/Dropbox/UBC/TReX_paper/figures/modelChannel.png -m5 -nofontswap
% 
% 
% %% SO how long is it going to take to get up the estuary
% cm=m_colmap('blues',7);cm([1 end-1:end],:)=[]; % Change this!
% 
% f1=figure('units','centimeters','outerposition',...
%     [0.01 0.01 20 15],'color','w');
% ax1=axes('position',[.1 .1 0.8 0.8]); hold on;
% lat_lim=[46 52];
% lon_lim=[-70 -56];
% m_proj('equidistant','lon',lon_lim,'lat',lat_lim);
% [CS,CH]=m_etopo2('contourf',0:-100:-600,'linestyle','none');
% caxis([-500 0]);
% colormap(cm);
% m_gshhs_l('patch',[0.9 0.9 0.9],'edgecolor','none');
% m_grid('linestyle','none','linewidth',0.1,'tickdir','out',...
%     'xaxisloc','top','yaxisloc','right','fontsize',6);
% s=m_scatter(dLon,dLat,200,'p','filled',...
%     'markerfacecolor',[0.929,0.694,0.125],'markeredgecolor','k');
% 
% lx=[-60.643,-61.7018,-62.5617,-63.7525,-64.8440,-64.8440,-66.0016,...
%     -66.8285,-67.4570,-68.1847];
% 
% 
% m_text(-58.0872-0.5,48.4215-0.4,'\itNewfoundland','fontsize',4);
% % m_text(-63.0546-0.3,49.4768,'\it{Anticosti Island}','fontsize',6);
% m_text(-66.4483-1.25,47.5129-.5,'\it{New Brunswick}','fontsize',4);
% % m_text(-63.5096,46.3486,'\it{P.E.I.}','fontsize',6);
% m_text(-66.2018,48.6351,'\it{Québec}','fontsize',4);
% m_text(-60,49.458724,{'Gulf of';'St. Lawrence'},'fontsize',4,...
%     'FontAngle', 'italic','HorizontalAlignment','center','color','w');
% m_text(-59.2638,47.1793,{'Cabot';'Strait'},'fontsize',4,...
%     'FontAngle', 'italic','HorizontalAlignment','center','color','w');
% 
% m_northarrow(-65.99,51.15,1,'type',3,'linewi',.5);
% 
% 
% %% 3D plot of Gulf with 275m slice
% 
% figure;
% surf(bathy.lon,bathy.lat,bathy.Z,'edgecolor','none');
% hold on
% caxis([-500 0]);
% colormap(m_colmap('blues'))
% o=ones(size(bathy.Z))*-275;
% mesh(bathy.lon,bathy.lat,o,'facealpha',0.2,...
%     'edgecolor','none','facecolor',rgb_x('pink'));
% 
% %% What % of area deeper than 200m is *NOT* the Laurentian channel
% load GoSL_Az.mat Az Z
% 
% LC.lon=[LC.lon -67.3251369175118,-67.8203816027231,-68.5137241620189,-68.9429362225353,...
%     -69.2730993460095,-69.6032624694837,-69.5702461571363,-69.2730993460095,-68.9759525348827,...
%     -68.5137241620189,-67.9854631644602,-67.5562511039437,-66.9949737940376,-66.4667127964789]; 
% 
% LC.lat=[LC.lat 49.5198555956679,49.3465703971119,49.1299638989170,48.9133574007220,...
%     48.6967509025271,48.3501805054152,48.0469314079422,48.0252707581227,48.1768953068592,...
%     48.3285198555957,48.5451263537906,48.6967509025271,48.8483754512635,49.0216606498195];
% LC.station_bound=alphaShape(LC.lon',LC.lat',1,'HoleThreshold',15);
%     
% Bdist=mode(m_lldist(bathy.lon(:),bathy.lat(:)))*1000; % m
% 
% idxN=find(~inShape(LC.station_bound,bathy.lon,bathy.lat) &...
%     bathy.lat>48 & bathy.lon>-63 & bathy.Z<=-250);
% Az_notLC=zeros(max(Z)+1,1);
% allZnotLC=bathy.Z(idxN)*-1;
% 
% idxLC=find(inShape(LC.station_bound,bathy.lon,bathy.lat)& bathy.Z<=-250);
% Az_LC=zeros(max(Z)+1,1);
% allZLC=bathy.Z(idxLC)*-1;
% 
% for i=1:length(Z)-1
%         Az_notLC(i+1)=sum(allZnotLC>i)*Bdist*Bdist; %m^2
%         Az_LC(i+1)=sum(allZLC>i)*Bdist*Bdist; %m^2
% end
% 
% volLC=sum(Az_LC(200:350))
% volNotLC=sum(Az_notLC(200:350))
% 
% volLC/volNotLC
% 
% figure;hold on;
% scatter(bathy.lon(idxLC),bathy.lat(idxLC));
% scatter(bathy.lon(idxN),bathy.lat(idxN));
% 
% % Area below 250m
% Az_LC(250)/(Az_LC(250)+Az_notLC(250))*100
% 
% %% junk code that might be useful
% % jbfill(centredX1(~isnan(ICmn1)),ICmn1(~isnan(ICmn1))-ICstd1(~isnan(ICmn1)),...
% %     flipud(ICmn1(~isnan(ICmn1))+ICstd1(~isnan(ICmn1))),...
% %     rgb_x('light blue'),rgb('blue'),1,0);
% % errorbar(centredX1,ICmn1,ICstd1,'markersize',5,'marker','o',...
% %     'markerfacecolor','b','markeredgecolor','k',...
% %     'color','b');
% 
% 
% % jbfill(centredX2(~isnan(ICmn2)),ICmn2(~isnan(ICmn2))-ICstd2(~isnan(ICmn2)),...
% %     flipud(ICmn2(~isnan(ICmn2))+ICstd2(~isnan(ICmn2))),...
% %     rgb_x('light red'),rgb('red'),1,0.2);
% % errorbar(centredX2,ICmn2,ICstd2,'markersize',5,'marker','o',...
% %     'markerfacecolor','r','markeredgecolor','k',...
% %     'color','r');
%    
% 
% % % find isopycnal injection sample
% % IP_TR1=NaN(length(idx1),1);IP_TRx1=NaN(length(idx1),1);
% % TR1=[];TRx1=[];
% % for i=1:length(idx1)
% %     [d,IPi]=min(abs(CTD1(idx1(i)).TR_sigma-27.26));
% %     TR1=[TR1;CTD1(idx1(i)).TR_SF5];
% %     TRx1=[TRx1;ones(length(CTD1(idx1(i)).TR_SF5),1)*stdist1(i)];
% %     if d<0.3
% %         IP_TR1(i)=CTD1(idx1(i)).TR_SF5(IPi);
% %         IP_TRx1(i)=stdist1(idx1(i));
% %     end
% % end
% % 
% % IP_TR2=NaN(length(idx2),1);IP_TRx2=NaN(length(idx2),1);
% % 
% % TR2=[];TRx2=[];
% % % find isopycnal injection sample
% % for i=1:length(idx2)
% %     [d,IPi]=min(abs(CTD2(idx2(i)).TR_sigma-27.26));
% %     TR2=[TR2;CTD2(idx2(i)).TR_SF5];
% %     TRx2=[TRx2;ones(length(CTD2(idx2(i)).TR_SF5),1)*stdist2(i)];
% %     if d<0.3
% %         IP_TR2(i)=CTD2(idx2(i)).TR_SF5(IPi);
% %         IP_TRx2(i)=stdist2(idx2(i));
% %     end
% % end
%     
%     
%     
%     
%     
%     
% 
% 
% 
% 
