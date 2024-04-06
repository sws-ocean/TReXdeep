%% Modelling analysis
addpath(genpath('/Users/samst/Dropbox/UBC/TReX_paper/'));
addpath(genpath('/Users/samst/Dropbox/UBC/TReX_deepII/'));
addpath(genpath('/Users/samst/Dropbox/UBC/GSW/'));
addpath(genpath('/Users/samst/Dropbox/UBC/m_map/'));
addpath(genpath('/Users/samst/Dropbox/UBC/Misc/'));

clear

%% load DATA
load('data20230707.mat', 'CTD1','CTD2','stlon1','stlon2','stlat1','stlat2')
MT.L1=load('/Users/samst/Library/CloudStorage/Dropbox/UBC/TReX_paper/code/letter/R1/conc_3d/conc_20220617_K10.mat');
MT.L2=load('/Users/samst/Library/CloudStorage/Dropbox/UBC/TReX_paper/code/letter/R1/conc_3d/conc_20221102_K10.mat');


% Load bathy and pull out station locs
bathy=load('GoSL_Zhires_mod.mat');
bathy.Z(bathy.Z==-9999)=NaN;
% Decimate bathy to speed things up
dec=10;
bathy.lon=bathy.lon(1:dec:end,1:dec:end);
bathy.lat=bathy.lat(1:dec:end,1:dec:end);
bathy.Z=bathy.Z(1:dec:end,1:dec:end);
model_z=griddata(bathy.lon,bathy.lat,bathy.Z,MT.L1.lon,MT.L1.lat);

%%
MT.L1.Chorz=sum(MT.L1.C,3,'omitnan'); %
MT.L2.Chorz=sum(MT.L2.C,3,'omitnan'); %

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
idx=find(inShape(LC.station_bound,MT.L1.lon,MT.L1.lat) & model_z<-200);

M1=sum(MT.L1.Chorz(find(inShape(LC.station_bound,MT.L1.lon,MT.L1.lat))));
M2=sum(MT.L2.Chorz(find(inShape(LC.station_bound,MT.L2.lon,MT.L2.lat))));

dLat=48.244; dLon=-59.952;

lat_lim=[46 51];
lon_lim=[-70 -57.5];

% m_proj('equidistant','lon',lon_lim,'lat',lat_lim);
% f1=figure('units','centimeters','outerposition',...
%     [0.01 0.01 20 20],'color','w');
% ax1=axes; hold on;%('position',[.0 .0 .1 1]); hold on;
% m_proj('equidistant','lon',lon_lim,'lat',lat_lim);
% [c,h]=m_contour(bathy.lon,bathy.lat,bathy.Z,0:-100:-500,'color',...
%     [0.95 0.95 0.95],'linewidth',0.5);
% [c,h]=m_contour(bathy.lon,bathy.lat,bathy.Z,[-279 -279],'color',...
%     'k','linewidth',0.5);
%
% clabel(c,h,'LabelSpacing',500,'color',[0.95 0.95 0.95],'fontsize',6);
% m_gshhs_i('patch',[0.9 0.9 0.9],'edgecolor','none');
% m_grid('linestyle','none','linewidth',0.1,'tickdir','out',...
%     'xaxisloc','top','yaxisloc','right','fontsize',6);
% s=m_scatter(dLon,dLat,200,'p','filled',...
%     'markerfacecolor',[0.929,0.694,0.125],'markeredgecolor','k');
%
% m_scatter(MT.L1.lon(idx),MT.L1.lat(idx),1,'x');
%
% [l,ll]=m_ll2xy(LC.lon,LC.lat);
% tmpshp=alphaShape(l',ll',0.01);
% LC.AS=plot(tmpshp,'facealpha',0.1,'edgecolor','none','facecolor',rgb_x('forest green'));

MT.L1.grids=load("/Users/samst/Library/CloudStorage/Dropbox/UBC/TReX_paper/code/letter/R1/ciopse_for_sam/modelGridsL1.mat");
MT.L2.grids=load("/Users/samst/Library/CloudStorage/Dropbox/UBC/TReX_paper/code/letter/R1/ciopse_for_sam/modelGridsL2.mat");

% station-wise analysis
MT.L1.stnIdxs=NaN(length(stlat1),1);

for i=1:length(MT.L1.stnIdxs)
    distances = sqrt((MT.L1.grids.lon(:) - stlon1(i)).^2 + (MT.L1.grids.lat(:) - stlat1(i)).^2);
    [~, MT.L1.stnIdxs(i,1)] = min(distances);
end

MT.L2.stnIdxs=NaN(length(stlat2),1);

for i=1:length(MT.L2.stnIdxs)
    distances = sqrt((MT.L2.grids.lon(:) - stlon2(i)).^2 + (MT.L2.grids.lat(:) - stlat2(i)).^2);
    [~, MT.L2.stnIdxs(i,1)] = min(distances);
end

%% lateral

flds=fieldnames(MT);

for i=1:2
    MT.(flds{i}).x_dist=gsw_distance([ones(size(MT.(flds{i}).lon(:)))*dLon MT.(flds{i}).lon(:)]...
        ,[ones(size(MT.(flds{i}).lat(:)))*dLat MT.(flds{i}).lat(:)])/1000;

    mDist=NaN(size(MT.(flds{i}).stnIdxs,1),1);
    mChannel=NaN(size(MT.(flds{i}).stnIdxs,1),1);

    for ii=1:length(MT.(flds{i}).stnIdxs)
        [I2(ii),J2(ii)] = ind2sub(size(MT.(flds{i}).lon),MT.(flds{i}).stnIdxs(ii));

        mChannel(ii,1)=sum(MT.(flds{i}).Chorz(I2(ii)-5:I2(ii)+5,...
            J2(ii)-5:J2(ii)+5),1:2,'omitnan');

        mDist(ii,1)=MT.(flds{i}).x_dist(MT.(flds{i}).stnIdxs(ii));
    end

    % % mChannel=MT.(flds{i}).Chorz(idx3);
    % mDist=MT.(flds{i}).x_dist(idx3);

    mDist(isnan(mChannel))=[];
    mChannel(isnan(mChannel))=[];

    % Make in-channel bin average
    interM=50;
    channelXM=0:interM/2:700;
    centredXM=interM/2:interM/2:max(channelXM)+(interM/2);
    count=0;
    MT.(flds{i}).ICmnM=NaN(length(channelXM),1);
    MT.(flds{i}).ICstdM=NaN(length(channelXM),1);

    for ii=channelXM
        count=count+1;
        MT.(flds{i}).ICmnM(count)=nanmean(mChannel(mDist>ii & mDist<ii+interM));
        MT.(flds{i}).ICstdM(count)=std(mChannel(mDist>ii & mDist<ii+interM),'omitnan');
    end

    MT.(flds{i}).mDist=mDist;
    MT.(flds{i}).mChannel=mChannel;
end

MT.L1.ICmnM(isnan(MT.L1.ICmnM))=0;
MT.L2.ICmnM(isnan(MT.L2.ICmnM))=0;

% can I model? yes!
M0=3e5; % g
X=0:1e3:725e3; % m
A=75:2:300; %m2/s
U=100/86400:15/86400:750/86400;% m/s
t=[(CTD1(1).mtime(1)-datenum(2021,10,24))*86400;...
    (CTD2(1).mtime(1)-datenum(2021,10,24))*86400]; % s
mm_sf5cf3 = 196.06;

ychannel=75e3; zspread=1; % m
Achannel=ychannel;%*zspread; % m2, area of the channel*depth range of tracer

c1horz=NaN(length(X),length(U)*length(A));
cost1=NaN(1,length(U)*length(A));

c2horz=NaN(length(X),length(U)*length(A));
cost2=NaN(1,length(U)*length(A));

idxs=NaN(2,length(U)*length(A));
count=0;

[tmpDistM1,idxM1]=sort(MT.L1.mDist);
tmpTRM=MT.L1.mChannel;tmpTRM1=tmpTRM(idxM1);

[tmpDistM2,idxM2]=sort(MT.L2.mDist);
tmpTRM=MT.L2.mChannel;tmpTRM2=tmpTRM(idxM2);

inL1=trapz(tmpDistM1,tmpTRM1);
inL2=trapz(tmpDistM2,tmpTRM2);

normL1=(tmpTRM1-min(tmpTRM1))./(max(tmpTRM1)-min(tmpTRM1));
normL2=(tmpTRM2-min(tmpTRM2))./(max(tmpTRM2)-min(tmpTRM2));

costM1=[]; costM2=[];

for i=1:length(A)
    for ii=1:length(U)

        count=count+1;

        c1horz(:,count)=M0./sqrt(4*pi*A(i)*t(1)).*exp( -(X-U(ii)*t(1)).^2./(4*A(i)*t(1)));%...
        % ./1e3.*1e15./mm_sf5cf3/Achannel.*1026; % fmol/m2
        c1horz(:,count)=(c1horz(:,count)-min(c1horz(:,count)))/(max(c1horz(:,count))-min(c1horz(:,count))).*max(MT.L1.ICmnM);

        c2horz(:,count)=M0./sqrt(4*pi*A(i)*t(2)).*exp( -(X-U(ii)*t(2)).^2./(4*A(i)*t(2)));%...
        % ./1e3.*1e15./mm_sf5cf3/Achannel.*1026; % fmol/m2
        c2horz(:,count)=(c2horz(:,count)-min(c2horz(:,count)))/(max(c2horz(:,count))-min(c2horz(:,count))).*max(MT.L2.ICmnM);

        % costM1(count)=trapz(tmpDistM1,(interp1(X,c1horz(:,count),tmpDistM1)-...
        %     tmpTRM1).^2);

        costM1(count)=trapz(centredXM*1000,(interp1(X,c1horz(:,count),centredXM*1000)-...
            MT.L1.ICmnM').^2);

        costM2(count)=trapz(centredXM*1000,(interp1(X,c2horz(:,count),centredXM*1000)-...
            MT.L2.ICmnM').^2);

        idxs(:,count)=[i;ii];
    end
end

[~,idxc1horz]=min(costM1); [~,idxc2horz]=min(costM2);

fprintf('\nLeg 1 1D model: k=%3.0f m^2 s^{-1}, U=%2.2f cm s^{-1}\n',...
    A(idxs(1,idxc1horz)),U(idxs(2,idxc1horz))*100);
fprintf('\nLeg 2 1D model: k=%3.0f m^2 s^{-1}, U=%2.2f cm s^{-1}\n',...
    A(idxs(1,idxc2horz)),U(idxs(2,idxc2horz))*100);

figure; hold on;
scatter(MT.L1.mDist,MT.L1.mChannel,10,'markerfacecolor','b',...
    'MarkerEdgeColor','none','MarkerFaceAlpha',0.05)
plot(centredXM,MT.L1.ICmnM,'b--','linewidth',1);
plot(X/1000,c1horz(:,idxc1horz),'b','linewidth',2)

scatter(MT.L2.mDist,MT.L2.mChannel,10,'markerfacecolor','r',...
    'MarkerEdgeColor','none','MarkerFaceAlpha',0.05)
plot(centredXM,MT.L2.ICmnM,'r--','linewidth',1);
plot(X/1000,c2horz(:,idxc2horz),'r','linewidth',2)
% ylim([0 500]);


%% errors
B1=bootrnd(size(stlon1,1),1000);
B2=bootrnd(size(stlon2,1),1000);

Aerr1=NaN(size(B1,2),1); Uerr1=NaN(size(B1,2),1);
Aerr2=NaN(size(B2,2),1); Uerr2=NaN(size(B2,2),1);

hw=waitbar(0,'Calculating...');
tic
for j=1:size(B1,2)

    idxsE=NaN(2,length(U)*length(A));
    count=0;

    cost1E=NaN(1,length(U)*length(A));
    cost2E=NaN(1,length(U)*length(A));

    c1E=NaN(length(X),length(U)*length(A));
    c2E=NaN(length(X),length(U)*length(A));

    for i=1:length(A)
        for ii=1:length(U)
            count=count+1;

            D1=MT.L1.mDist(B1(:,j));
            T1=MT.L1.mChannel(B1(:,j));
            [D1,idxD1]=sort(D1);
            T1=T1(idxD1);

            ICmnM=[];
            for iii=1:length(channelXM)
                ICmnM(iii)=nanmean(T1(D1>channelXM(iii) & D1<channelXM(iii)+interM));
            end

            c1E(:,count)=M0./sqrt(4*pi*A(i)*t(1)).*exp( -(X-U(ii)*t(1)).^2./(4*A(i)*t(1)));%...
            c1E(:,count)=(c1E(:,count)-min(c1E(:,count)))/...
                (max(c1E(:,count))-min(c1E(:,count))).*max(ICmnM);

            cost1E(count)=trapz(D1*1000,(interp1(X,c1E(:,count),D1*1000)-...
                T1).^2);

            D2=MT.L2.mDist(B2(:,j));
            T2=MT.L2.mChannel(B2(:,j));
            [D2,idxD2]=sort(D2);
            T2=T2(idxD2);

            ICmnM=[];
            for iii=1:length(channelXM)
                ICmnM(iii)=nanmean(T2(D2>channelXM(iii) & D2<channelXM(iii)+interM));
            end

            c2E(:,count)=M0./sqrt(4*pi*A(i)*t(2)).*exp( -(X-U(ii)*t(2)).^2./(4*A(i)*t(2)));%...
            c2E(:,count)=(c2E(:,count)-min(c2E(:,count)))/...
                (max(c2E(:,count))-min(c2E(:,count))).*max(ICmnM);

            cost2E(count)=trapz(D2*1000,(interp1(X,c2E(:,count),D2*1000)-...
                T2).^2);


            idxsE(:,count)=[i;ii];
        end
    end
    [~,idxC1E]=min(cost1E); [~,idxC2E]=min(cost2E);

    Aerr1(j,1)=A(idxsE(1,idxC1E));
    Uerr1(j,1)=U(idxsE(2,idxC1E));
    Aerr2(j,1)=A(idxsE(1,idxC2E));
    Uerr2(j,1)=U(idxsE(2,idxC2E));

    waitbar(j/size(B1,2),hw);
end
toc
delete(hw)

% From 1000 subsets
fprintf('\nLeg 1 (2 STD): k=%3.0f m^2 s^{-1}, U=%4.3f cm s^{-1}\n',...
    std(Aerr1)*2,std(Uerr1)*2*100);
fprintf('\nLeg 2 (2 STD): k=%3.0f m^2 s^{-1}, U=%4.3f cm s^{-1}\n',...
    std(Aerr2)*2,std(Uerr2)*2*100);

%% Vertical

MT.L1.grids=load("/Users/samst/Library/CloudStorage/Dropbox/UBC/TReX_paper/code/letter/R1/ciopse_for_sam/modelGridsL1.mat");
MT.L2.grids=load("/Users/samst/Library/CloudStorage/Dropbox/UBC/TReX_paper/code/letter/R1/ciopse_for_sam/modelGridsL2.mat");


Lmsk=MT.L1.lon(:,1)>-69 & MT.L1.lon(:,1)<-56;
LLmsk=MT.L1.lat(1,:)>46 & MT.L1.lat(1,:)<51;

MT.L1.Csub=MT.L1.C(find(Lmsk,1,'first'):find(Lmsk,1,'last'),...
    find(LLmsk,1,'first'):find(LLmsk,1,'last'),:);

MT.L2.Csub=MT.L2.C(find(Lmsk,1,'first'):find(Lmsk,1,'last'),...
    find(LLmsk,1,'first'):find(LLmsk,1,'last'),:);

% Take mean sigma from CS casts
load('data20230707.mat', 'CScast1','CScast2','sigCo1','sigCo2');

stSigM1=NaN(450,length(stlat1));
stSigM2=NaN(450,length(stlat1));

% station-wise analysis
stnIdxs1=NaN(length(stlat1),1);

for i=1:length(stnIdxs1)
    distances = sqrt((MT.L1.grids.lon(:) - stlon1(i)).^2 + (MT.L1.grids.lat(:) - stlat1(i)).^2);
    [~, stnIdxs1(i,1)] = min(distances);

    [I1(i),J1(i)] = ind2sub(size(MT.L1.grids.lon),stnIdxs1(i));
    stSigM1(:,i)=interp1(MT.L1.grids.modelZ,...
        squeeze(MT.L1.grids.modelSigma(I1(i),J1(i),:)),1:450);
end

stnIdxs2=NaN(length(stlat2),1);

for i=1:length(stnIdxs2)
    distances = sqrt((MT.L2.grids.lon(:) - stlon2(i)).^2 + (MT.L2.grids.lat(:) - stlat2(i)).^2);
    [~, stnIdxs2(i,1)] = min(distances);

    [I2(i),J2(i)] = ind2sub(size(MT.L1.grids.lon),stnIdxs2(i));
    stSigM2(:,i)=interp1(MT.L2.grids.modelZ,...
        squeeze(MT.L2.grids.modelSigma(I2(i),J2(i),:)),1:450);
end

figure;hold on;
scatter(MT.L1.grids.lon(:),MT.L1.grids.lat(:),2)
scatter(stlon1,stlat1,10,'x');
scatter(MT.L1.grids.lon(stnIdxs1),MT.L1.grids.lat(stnIdxs1),40,'o');
scatter(stlon2,stlat2,10,'x');
scatter(MT.L1.grids.lon(stnIdxs2),MT.L1.grids.lat(stnIdxs2),40,'^');

sigCoM1=mean(stSigM1(:,CScast1),2,'omitnan');
sigCoM2=mean(stSigM2(:,CScast2),2,'omitnan');

figure; hold on;
plot(sigCo1,1:450,'b');
plot(sigCoM1,1:450,'b--');

plot(sigCo2,1:450,'r');
plot(sigCoM2,1:450,'r--');
axis ij
grid on

zStar_grid=[1:450]';
zHat_grid=[-300:1:200]';

CijM1=NaN(length(zHat_grid),length(stlat1));
CijM2=NaN(length(zHat_grid),length(stlat2));
zStarM1=NaN(length(zStar_grid),length(stlat1));
zStarM2=NaN(length(zStar_grid),length(stlat2));
zHatM1=NaN(length(zStar_grid),length(stlat1));
zStarM2=NaN(length(zStar_grid),length(stlat2));

Gspacing=median(gsw_distance(MT.L1.grids.lon(:),MT.L1.grids.lat(:)));

% Leg 1 station wise
for i=1:length(stlon1)

    % subset virtual tracer cast within 5 km
    tmpTR=squeeze(sum(MT.L1.Csub(I1(i)-3:I1(i)+3,...
        J1(i)-3:J1(i)+3,:),1:2));
    tmpSigma=stSigM1(:,i);

    % interpolate these onto CS mean sigma/zHat_grid relationship
    zStarM1(:,i)=interp1(sigCoM1,zStar_grid,tmpSigma);
    zHatM1(:,i)=zStarM1(:,i)-279;
    msk=~isnan(zHatM1(:,i));

    CijM1(:,i)=interp1(zHatM1(msk,i),tmpTR(msk)',zHat_grid);
    CijM1(isnan(CijM1(:,i)),i)=0;
end

% Leg 2 station wise
for i=1:length(stlon2)

    % subset virtual tracer cast within 5 km
    tmpTR=squeeze(sum(MT.L2.Csub(I2(i)-3:I2(i)+3,...
        J2(i)-3:J2(i)+3,:),1:2));
    tmpSigma=stSigM2(:,i);

    % interpolate these onto CS mean sigma/zHat_grid relationship
    zStarM2(:,i)=interp1(sigCoM2,zStar_grid,tmpSigma);
    zHatM2(:,i)=zStarM2(:,i)-279;
    msk=~isnan(zHatM2(:,i));

    CijM2(:,i)=interp1(zHatM2(msk,i)+(rand(size(zHatM2(msk,i)))/100),tmpTR(msk)',zHat_grid);
    CijM2(isnan(CijM2(:,i)),i)=0;
end

figure;hold on;
plot(CijM1,zHat_grid,'color',rgb_x('cornflower'),'LineWidth',0.1);
plot(CijM2,zHat_grid,'color',rgb_x('dusty pink'),'LineWidth',0.1);
plot(mean(CijM1,2),zHat_grid,'b','linewidth',2);
plot(mean(CijM2,2),zHat_grid,'r','linewidth',2);
ylim([-300 300]);xlim([0 5]);
axis ij;

%% Advection-diffusion model
M0=[3e5 3e5]; % n
A0=[2.1e11 2.1e11] ;% m^2
L=200; % m
kappa=[0.5e-6:2e-6:1e-3].*86400; % m^2/day
U=-150/365/86400:1/365/86400:25/365/86400; % m/s

t=[(CTD1(1).mtime(1)-datenum(2021,10,24)) ...
    (CTD2(1).mtime(1)-datenum(2021,10,24))];  % day

idxs=NaN(2,length(U)*length(kappa));

% run through different kappas and find the best fit
cost1=NaN(length(kappa)*length(U),1);
cost2=NaN(length(kappa)*length(U),1);

c1=NaN(length(zHat_grid),length(kappa)*length(U));
c2=NaN(length(zHat_grid),length(kappa)*length(U));

CM1=mean(CijM1,2);
CM2=mean(CijM2,2);

riM1=CM1/sum(CM1);
riM2=CM2/sum(CM2);

count=0;
hw=waitbar(0,'Calculating...');
for i=1:length(kappa)
    for ii=1:length(U)
        count=count+1;

        w=U(ii)*86400;

        W=w+kappa(i)/L; % m/s

        c1(:,count)=M0(1)/A0(1)/sqrt(4*pi*kappa(i)*t(1))*...
            exp(-(zHat_grid-W*t(1)).^2./(4*kappa(i).*t(1)) + w*t(1)/L);
        % c1(:,count)=(c1(:,count)-min(c1(:,count)))/(max(c1(:,count))-min(c1(:,count))).*(max(smooth(mean(CijM1,2),10)).*1);
        c1(:,count)=c1(:,count)/sum(c1(:,count));


        c2(:,count)=M0(2)/A0(2)/sqrt(4*pi*kappa(i)*t(2))*...
            exp(-(zHat_grid-W*t(2)).^2./(4*kappa(i).*t(2)) + w*t(2)/L);
        % c2(:,count)=(c2(:,count)-min(c2(:,count)))/(max(c2(:,count))-min(c2(:,count))).*(max(smooth(mean(CijM2,2),10)).*1);
        c2(:,count)=c2(:,count)/sum(c2(:,count));

        cost1(count)=trapz(zHat_grid,(c1(:,count)-riM1).^2);
        cost2(count)=trapz(zHat_grid,(c2(:,count)-riM2).^2);

        idxs(:,count)=[i;ii];

    end
    waitbar(i/length(kappa),hw);
end

% find minimum cost
[~,Cidx1]=min(cost1);
[~,Cidx2]=min(cost2);

clc
fprintf('Leg 1 1D model: k=%2.1e m^2 s^{-1}, w=%2.1f m yr^{-1}\n',...
    kappa(idxs(1,Cidx1))/86400,U(idxs(2,Cidx1))*365*86400);
fprintf('Leg 2 1D model: k=%2.1e m^2 s^{-1}, w=%2.1f m yr^{-1}\n',...
    kappa(idxs(1,Cidx2))/86400,U(idxs(2,Cidx2))*365*86400);

% figure;hold on;
% plot(CijM1,zHat_grid,'color',rgb_x('cornflower'),'LineWidth',0.1);
% plot(CijM2,zHat_grid,'color',rgb_x('dusty pink'),'LineWidth',0.1);
% plot(mean(CijM1,2),zHat_grid,'b--','linewidth',2);
% plot(mean(CijM2,2),zHat_grid,'r--','linewidth',2);
% plot(c1(:,Cidx1),zHat_grid,'b','linewidth',2)
% plot(c2(:,Cidx2),zHat_grid,'r','linewidth',2)
% axis ij;
% ylim([-300 300]);xlim([0 5]);
% axis ij;

figure;hold on;
plot(riM1,zHat_grid,'b--','linewidth',2);
plot(riM2,zHat_grid,'r--','linewidth',2);
plot(c1(:,Cidx1),zHat_grid,'b','linewidth',2)
plot(c2(:,Cidx2),zHat_grid,'r','linewidth',2)
axis ij;
ylim([-300 300]); grid on;
axis ij;


%% Boundary vs interior
load('data20230707.mat', 'stBath1','stBath2');

tic
% Leg 1
rM1In=mean(CijM1(:,stBath1<-279),2)./sum(mean(CijM1(:,stBath1<-279),2));
cost1InDiff=c1-rM1In;
cost1In=trapz(zHat_grid,cost1InDiff.^2);
[~,Cidx1In]=min(cost1In);

rM1Bn=mean(CijM1(1:300,stBath1>-279),2)./sum(mean(CijM1(1:300,stBath1>-279),2));
cost1BnDiff=c1(1:300,:)-rM1Bn;
cost1Bn=trapz(zHat_grid(1:300),cost1BnDiff.^2);
[~,Cidx1Bn]=min(cost1Bn);

% Leg 2
rM2In=mean(CijM2(:,stBath2<-279),2)./sum(mean(CijM2(:,stBath2<-279),2));
cost2InDiff=c2-rM2In;
cost2In=trapz(zHat_grid,cost2InDiff.^2);
[~,Cidx2In]=min(cost2In);

rM2Bn=mean(CijM2(1:300,stBath2>-279),2)./sum(mean(CijM2(1:300,stBath2>-279),2));
cost2BnDiff=c2(1:300,:)-rM2Bn;
cost2Bn=trapz(zHat_grid(1:300),cost2BnDiff.^2);
[~,Cidx2Bn]=min(cost2Bn);
toc

fprintf('Leg 1 INTERIOR: k=%2.1e m^2 s^{-1}, w=%2.1f m yr^{-1}\n',...
    kappa(idxs(1,Cidx1In))/86400,U(idxs(2,Cidx1In))*365*86400);
fprintf('Leg 1 BOUNDARY: k=%2.1e m^2 s^{-1}, w=%2.1f m yr^{-1}\n',...
    kappa(idxs(1,Cidx1Bn))/86400,U(idxs(2,Cidx1Bn))*365*86400);

fprintf('Leg 2 INTERIOR: k=%2.1e m^2 s^{-1}, w=%2.1f m yr^{-1}\n',...
    kappa(idxs(1,Cidx2In))/86400,U(idxs(2,Cidx2In))*365*86400);
fprintf('Leg 2 BOUNDARY: k=%2.1e m^2 s^{-1}, w=%2.1f m yr^{-1}\n',...
    kappa(idxs(1,Cidx2Bn))/86400,U(idxs(2,Cidx2Bn))*365*86400);

figure;
subplot(1,2,1);hold on;
plot(rM1In,zHat_grid,'b--','linewidth',1);
plot(rM2In,zHat_grid,'r--','linewidth',1);
plot(c1(:,Cidx1In),zHat_grid,'color',rgb_x('cornflower'),'linewidth',2)
plot(c2(:,Cidx2In),zHat_grid,'color',rgb_x('dusty pink'),'linewidth',2)
axis ij;
ylim([-300 300]); title('Interior');

subplot(1,2,2);hold on;
plot(rM1Bn(1:300),zHat_grid(1:300),'b:','linewidth',1);
plot(rM2Bn(1:300),zHat_grid(1:300),'r:','linewidth',1);
plot(c1(:,Cidx1Bn),zHat_grid,'color',rgb_x('cornflower'),'linewidth',2)
plot(c2(:,Cidx2Bn),zHat_grid,'color',rgb_x('dusty pink'),'linewidth',2)
axis ij;
ylim([-300 300]); title('Boundary');


% The vertical spreading analysis shows enhanced vertical spreading an
% order of magnitude larger than Observed, suggesting that vertical
% mixing is too strong. When subsetting the data, we see enhanced shoaling
% at the boundaries, but a roughly consistent kappa_z throughout the basin,
% suggesting that the asymettircal nature of the mixing is not well
% simulated in the model, or, more likely, that the boundary effects of
% mixing are strongly incorporated into the basin-wide spreading during
% both Leg 1 and 2.

%% entire basin

Ci1=squeeze(mean(MT.L1.Csub,1:2,'omitnan'));Ci1(1)=NaN;
ri1=Ci1./nansum(Ci1);
ri1=interp1([1:500]-279,ri1,zHat_grid);
ri1(isnan(ri1))=0;
Ci2=squeeze(mean(MT.L2.Csub,1:2,'omitnan'));Ci2(1)=NaN;
ri2=Ci1./nansum(Ci2);
ri2=interp1([1:500]-275,ri2,zHat_grid);
ri2(isnan(ri2))=0;

% Leg 1
cost1Diff=c1-ri1;
cost1All=trapz(zHat_grid,cost1Diff.^2);
[~,Cidx1All]=min(cost1All);

% Leg 2
cost2Diff=c2-ri2;
cost2All=trapz(zHat_grid,cost2Diff.^2);
[~,Cidx2All]=min(cost2All);

fprintf('Leg 1 ALL: k=%2.1e m^2 s^{-1}, w=%2.1f m yr^{-1}\n',...
    kappa(idxs(1,Cidx1All))/86400,U(idxs(2,Cidx1All))*365*86400);

fprintf('Leg 2 ALL: k=%2.1e m^2 s^{-1}, w=%2.1f m yr^{-1}\n',...
    kappa(idxs(1,Cidx2All))/86400,U(idxs(2,Cidx2All))*365*86400);

figure;hold on;
plot(ri1,zHat_grid,'b--','linewidth',2);
plot(ri2,zHat_grid,'r--','linewidth',2);
plot(c1(:,Cidx1All),zHat_grid,'b','linewidth',2)
plot(c2(:,Cidx2All),zHat_grid,'r','linewidth',2)
axis ij;
ylim([-300 300]); grid on;
axis ij;


%% Mapping
lat_lim=[46.75-2 52];
lon_lim=[-69.75 -57.5+1];
col=[255 214 140]/255; % YELLOW!
f1=figure('units','centimeters','outerposition',...
    [0.01 0.01 21 25],'color','w');

%%%%%% maps %%%%%%
%%% L1 %%
ax1=axes('position',[0.1 0.6 0.4 0.4]);
hold on;
m_proj('equidistant','lon',lon_lim,'lat',lat_lim);
m_gshhs_i('patch',col,'edgecolor',[0.3 0.3 0.3],'linewidth',0.3);
m_grid('linestyle','none','linewidth',0.1,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','left','fontsize',8);
lev=1:5:120;
[c,h]=m_contourf(MT.L1.lon,MT.L1.lat,MT.L1.Chorz,lev,...
    'linestyle','none');
tmp=seminfhaxby;tmp2=linspecer;
colormap([[0.9 0.9 0.9];tmp2]);
clim([1 120]);

m_contour(bathy.lon,bathy.lat,bathy.Z,[-275 -275],...
    'color',rgb_x('dark grey'),'linewidth',0.66,'linestyle','-');
hh=m_contfbar(gca,0.175,[0.2 0.4],c,h,'endpiece','yes','xaxisloc','top',...
    'axfrac',0.02,'edgecolor','none');

ylabel(hh,{'\it{P_j}'});
s=m_scatter(dLon,dLat,200,'p','filled',...
    'markerfacecolor',rgb_x('light violet'),'markeredgecolor','k');
m_text(-69.3,45,'a)','fontsize',8);
title('Leg 1 (8 months)');

m_ruler([0.5 0.8],0.9,'fontsize',8,'linewidth',0.75);

%%%%%% L2 %%%%%%
ax1=axes('position',[0.525 0.6 0.4 0.4]);
hold on;
m_proj('equidistant','lon',lon_lim,'lat',lat_lim);
m_gshhs_i('patch',col,'edgecolor',[0.3 0.3 0.3],'linewidth',0.3);
m_grid('linestyle','none','linewidth',0.1,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','left','fontsize',8,'yticklabels',[]);
[c,h]=m_contourf(MT.L2.lon,MT.L2.lat,MT.L2.Chorz,lev,...
    'linestyle','none');
colormap([[0.9 0.9 0.9];tmp2]);
clim([1 120]);
m_contour(bathy.lon,bathy.lat,bathy.Z,[-275 -275],...
    'color',rgb_x('dark grey'),'linewidth',0.66,'linestyle','-');

s=m_scatter(dLon,dLat,200,'p','filled',...
    'markerfacecolor',rgb_x('light violet'),'markeredgecolor','k');

m_text(-69.3,45,'b)','fontsize',8);
title('Leg 2 (12 months)');

%%% Lateral %%%
ax3=axes('position',[0.1 0.5 0.325 0.1]); hold on;
title('Simulated');
plot(centredXM,MT.L1.ICmnM,'b--','linewidth',0.5);
scatter(centredXM,MT.L1.ICmnM,10,'marker','s','markerfacecolor','b',...
    'markeredgecolor','none','markerfacealpha',0.5);
plot(X/1000,c1horz(:,idxc1horz),'-','color',rgb_x('cornflower'),'linewidth',2)

plot(centredXM,MT.L2.ICmnM,'r--','linewidth',0.5);
scatter(centredXM,MT.L2.ICmnM,10,'marker','s','markerfacecolor','r',...
    'markeredgecolor','none','markerfacealpha',0.5);
plot(X/1000,c2horz(:,idxc2horz),'-','color',rgb_x('dusty pink'),'linewidth',2)
xticklabels([]);
yticks(0:1000:4000);
ylabel('\it{P_j}');
axis tight; grid on;
xlim([0 525]); xticks(0:100:500);
text(15,3,'c)','fontsize',8);

lp(1)=plot(NaN,NaN,'k--','marker','s','markerfacecolor','k',...
    'markeredgecolor','none');
lp(2)=plot(NaN,NaN,'k-','linewidth',1.5);

ax4=axes('position',[0.1 0.4 0.325 0.075]); hold on;
horzObs=load('/Users/samst/Library/CloudStorage/Dropbox/UBC/TReX_paper/code/letter/R1/channel1D_R1.mat');
title('Observed');
plot(horzObs.centredX1,horzObs.ICmn1,'b--','linewidth',0.5);
scatter(horzObs.centredX1,horzObs.ICmn1,10,'marker','s','markerfacecolor','b',...
    'markeredgecolor','none','markerfacealpha',0.5);
plot(horzObs.centredX2,horzObs.ICmn2,'r--','linewidth',0.5);
scatter(horzObs.centredX2,horzObs.ICmn2,10,'marker','s','markerfacecolor','r',...
    'markeredgecolor','none','markerfacealpha',0.5);

plot(horzObs.X/1e3,horzObs.cf1,'-','color',rgb_x('cornflower'),'linewidth',2);
plot(horzObs.X/1e3,horzObs.cf2,'-','color',rgb_x('dusty pink'),'linewidth',2);
ylabel('\it{I_j} \rm (fmol m^{-2})');
% ylim([0 max(horzObs.cf1)]);
xlim([0 525]); grid on;
xlabel('\it{x} \rm (km)');
text(15,0.3e4,'d)','fontsize',8);


% vertical
ax5=axes('position',[0.475 0.4 0.1 0.2]); hold on;
plot(riM1,zHat_grid,'b--','linewidth',0.5);
plot(riM2,zHat_grid,'r--','linewidth',0.5);
plot(c1(:,Cidx1),zHat_grid,'color',rgb_x('cornflower'),'linewidth',2)
plot(c2(:,Cidx2),zHat_grid,'color',rgb_x('dusty pink'),'linewidth',2)
axis ij;
ylim([-300 200]); grid on;
ylabel('$$\hat{z}\, (m)$$', 'Interpreter', 'latex');
xlabel('\it{r}');
title('Simulated');
text(0.008,-275,'e)');
vertObs=load('/Users/samst/Library/CloudStorage/Dropbox/UBC/TReX_paper/code/letter/R1/vertObs.mat');
xtickangle(0)

ax6=axes('position',[0.59 0.4 0.06 0.2]); hold on;
plot(vertObs.ri1,vertObs.zHat_grid,'b--','linewidth',0.5);
plot(vertObs.ri2,vertObs.zHat_grid,'r--','linewidth',0.5);
plot([0 0],[-300 -200],'color',rgb_x('cornflower'),'linewidth',2)
plot([0 0],[-300 -200],'color',rgb_x('dusty pink'),'linewidth',2)
plot(vertObs.cf1,vertObs.zHat_grid,'color',rgb_x('cornflower'),'linewidth',2)
plot(vertObs.cf2,vertObs.zHat_grid,'color',rgb_x('dusty pink'),'linewidth',1)
axis ij;
ylim([-300 200]); grid on;
yticklabels([]);
% xlabel('\it{r}');
title('Observed');
text(0.018,-275,'f)');
xtickangle(0)

% asymmetrical
ax7=axes('position',[0.475+0.25 0.4 0.1 0.2]); hold on;
plot(rM1In,zHat_grid,'--','color',rgb_x('cerulean'),'linewidth',0.5);
plot(rM2In,zHat_grid,'--','color',rgb_x('cerulean'),'linewidth',0.5);

plot(rM1Bn(1:300),zHat_grid(1:300),'--','color',...
    rgb_x('algae green'),'linewidth',0.5);
plot(rM2Bn(1:300),zHat_grid(1:300),'--','color',...
    rgb_x('algae green'),'linewidth',0.5);

l=plot(c1(:,Cidx1Bn),zHat_grid,'color',rgb_x('cornflower'),'linewidth',2);
l.Color=[rgb_x('algae green').*0.8 0.7];

l=plot(c2(:,Cidx2Bn),zHat_grid,'color',rgb_x('cornflower'),'linewidth',2);
l.Color=[rgb_x('algae green').*0.8 0.7];

l=plot(c1(:,Cidx1In),zHat_grid,'color',rgb_x('cornflower'),'linewidth',2);
l.Color=[rgb_x('cerulean') 0.8];

l=plot(c2(:,Cidx2In),zHat_grid,'color',rgb_x('cornflower'),'linewidth',2);
l.Color=[rgb_x('cerulean') 0.8];

axis ij;
ylim([-300 200]); grid on;
ylabel('$$\hat{z}\, (m)$$', 'Interpreter', 'latex');
xlabel('\it{r}');
title('Simulated');
text(0.018,-275,'g)');
xtickangle(0)

% l1=line([0.002 0.007],[110;110],'linewidth',2);
% l1.Color=[rgb_x('cerulean') 0.5];
% l2=plot([0.002 0.007],[130 130],'linewidth',2);
% l2.Color=[rgb_x('algae green') 0.5];

ax8=axes('position',[0.59+0.26 0.4 0.06 0.2]); hold on;
plot(vertObs.deep_ri1,vertObs.zHat_grid,'--','color',rgb_x('cerulean'),'linewidth',0.5);
plot(vertObs.deep_ri2,vertObs.zHat_grid,'--','color',rgb_x('cerulean'),'linewidth',0.5);

plot(vertObs.shallow_ri1,vertObs.zHat_grid,'--','color',...
    rgb_x('algae green'),'linewidth',0.5);
plot(vertObs.shallow_ri2,vertObs.zHat_grid,'--','color',...
    rgb_x('algae green'),'linewidth',0.5);

l=plot(vertObs.cfBn1,vertObs.zHat_grid,'color',rgb_x('cornflower'),'linewidth',2);
l.Color=[rgb_x('algae green').*0.8 0.7];

l=plot(vertObs.cfBn2,vertObs.zHat_grid,'color',rgb_x('cornflower'),'linewidth',2);
l.Color=[rgb_x('algae green').*0.8 0.7];

l=plot(vertObs.cfIn1,vertObs.zHat_grid,'color',rgb_x('cornflower'),'linewidth',2);
l.Color=[rgb_x('cerulean') 0.8];

l=plot(vertObs.cfIn1,vertObs.zHat_grid,'color',rgb_x('cornflower'),'linewidth',2);
l.Color=[rgb_x('cerulean') 0.8];

plot([0 0],[-300 -200],'color',rgb_x('cerulean'),'linewidth',2)

axis ij;
ylim([-300 200]); grid on;
title('Observed');
yticklabels([]);
text(0.018,-275,'h)');
xtickangle(0)

set(findall(gcf,'-property','fontsize'),'fontsize',8);

axes(ax7)
text(0.002,145,'Interior    (\it{h} \rm{>279 m)}','fontsize',6,'color',...
    rgb_x('cerulean'));
text(0.002,170,'Boundary (\it{h} \rm{<279 m)}','fontsize',6,'color',...
    rgb_x('algae green'));

axes(ax3)
legend(lp,'Window mean','1D model','location','northeast',...
    'fontsize',6);

axes(ax4)
text(410,2.5e4,'Leg 1','fontsize',6,'color',rgb_x('cornflower'));
text(410,2.1e4,'Leg 2','fontsize',6,'color',rgb_x('deep pink'));

set(findall(gcf,'-property','fontname'),'fontname','Latin Modern Roman');

%%
export_fig /Users/samst/Dropbox/UBC/TReX_paper/figures/R1/model_comp.pdf -dpdf -nofontswap