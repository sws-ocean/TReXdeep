%% Vertical mixing analysis
addpath(genpath('/Users/samst/Dropbox/UBC/TReX_deepII/'));
addpath(genpath('/Users/samst/Dropbox/UBC/TReX_paper/'));
addpath(genpath('/Users/samst/Dropbox/UBC/GSW/'));
addpath(genpath('/Users/samst/Dropbox/UBC/m_map/'));
addpath(genpath('/Users/samst/Dropbox/UBC/Misc/'));

clear
CTD1=load('TReX_DeepII_ctd.mat');CTD1=CTD1.CTD;
CTD2=load('TReX_DeepIV_ctd.mat');CTD2=CTD2.CTD;

% Load bathy and pull out station locs
bathy=load('GoSL_Zhires.mat');
bathy.Z(bathy.Z==-9999)=NaN;

% Decimate bathy to speed things up
dec=25;
bathy.lon=bathy.lon(1:dec:end,1:dec:end);
bathy.lat=bathy.lat(1:dec:end,1:dec:end);
bathy.Z=bathy.Z(1:dec:end,1:dec:end);

% Define MM of SF5
mm_sf5cf3 = 196.06;

% Define injection location
dLat=48.244; dLon=-59.952;

% Clean data and create data fields
CTD1(6).TR_SF5(1)=0; %remove spurious data
stlon1=NaN(length(CTD1),1);
stlat1=NaN(length(CTD1),1);
stz1=NaN(length(CTD1),1);
msk=zeros(length(CTD1),1);
for i=1:length(CTD1)
    CTD1(i).TR_SF5=CTD1(i).TR_SF5*1000; % Convert from pmol to fmol
    stlon1(i)=CTD1(i).longitude(1);
    stlat1(i)=CTD1(i).latitude(1);
    stz1(i)=max(CTD1(i).depSM);
    CTD1(i).x_dist=gsw_distance([dLon stlon1(i)],[dLat stlat1(i)])/1000;
    stBath1(i)=griddata(bathy.lon,bathy.lat,bathy.Z,...
        stlon1(i),stlat1(i),'nearest');
    % Remove stations with no SF5 data
    if isempty(CTD1(i).TR_SF5)
        msk(i)=1;
    end
end
CTD1(logical(msk))=[];
stlon1(logical(msk))=[];
stlat1(logical(msk))=[];
stz1(logical(msk))=[];
stBath1(logical(msk))=[];

stlon2=NaN(length(CTD2),1);
stlat2=NaN(length(CTD2),1);
stz2=NaN(length(CTD2),1);
msk=zeros(length(CTD2),1);
for i=1:length(CTD2)
    CTD2(i).TR_SF5=CTD2(i).TR_SF5*0.83;
    stlon2(i)=CTD2(i).longitude(1);
    stlat2(i)=CTD2(i).latitude(1);
    stz2(i)=max(CTD2(i).depSM);
    CTD2(i).x_dist=gsw_distance([dLon stlon2(i)],[dLat stlat2(i)])/1000;
    stBath2(i)=griddata(bathy.lon,bathy.lat,bathy.Z,...
        stlon2(i),stlat2(i));
    % Remove stations with no SF5 data
    if isempty(CTD2(i).TR_SF5)
        msk(i)=1;
    end
end
CTD2(logical(msk))=[];
stlon2(logical(msk))=[];
stlat2(logical(msk))=[];
stz2(logical(msk))=[];
stBath2(logical(msk))=[];

clearvars bathy

%% Start by interpolating tracer concs. onto isopycnal coordinates from CS

% Define CS shape
CS.lon=[-60.6601,-61.4159,-60.7996,-60.1602,-59.3463,-59.3928,-59.6253,-57.5674,...
    -58.7998,-58.0441,-59.8579];
CS.lat=[48.7870,48.0859,47.8163,47.2461,47.5851,48.1629,48.3401,46.5913,...
    46.0289,47.2153,46.6529];
CS.station_bound=alphaShape(CS.lon',CS.lat',1,'HoleThreshold',15);

% Find casts in CS
CScast1=find(inShape(CS.station_bound,stlon1,stlat1) & stz1>400);
CScast2=find(inShape(CS.station_bound,stlon2,stlat2) & stz2>400);

zStar_grid=[1:450]';
sigma_grid1=NaN(length(zStar_grid),length(CScast1));
sigma_grid2=NaN(length(zStar_grid),length(CScast2));

% Mean CS sigma for isopycnal coordinates
for i=1:length(CScast1)
    sigma_grid1(:,i)=interp1(CTD1(CScast1(i)).depSM,CTD1(CScast1(i)).sigma,zStar_grid);
end
sigCo1=mean(sigma_grid1,2,'omitnan');
[~,sigIdx1]=min(abs(sigCo1-27.26));

for i=1:length(CScast2)
    sigma_grid2(:,i)=interp1(CTD2(CScast2(i)).depSM,CTD2(CScast2(i)).sigma,zStar_grid);
end
sigCo2=mean(sigma_grid2,2,'omitnan');
[~,sigIdx2]=min(abs(sigCo2-27.26));

zHat_grid=[-200:1:200]';
lnes=lines(200);

G_or_L='linear'; % Gaussian or linear interpolation

for i=1:length(CTD1)
    
    % Need unique data for intepolation... add tiny variability
    tmpTR=CTD1(i).TR_SF5;
    tmpSigma=CTD1(i).TR_sigma+(rand(size(CTD1(i).TR_z))/100);
    
    % interpolate these onto CS mean sigma/zHat_grid relationship
    CTD1(i).zStar=interp1(sigCo1,zStar_grid,tmpSigma);
    
    % Change coordinates to distance from injection sigma
    CTD1(i).zHat=CTD1(i).zStar-zStar_grid(sigIdx1);
    
    % Remove any data outside of coordinate grid
    tmpTR(isnan(CTD1(i).zStar))=[];
    CTD1(i).zStar(isnan(CTD1(i).zStar))=[];
    CTD1(i).TR_SF5(isnan(CTD1(i).zHat))=[];
    CTD1(i).TR_z(isnan(CTD1(i).zHat))=[];
    CTD1(i).TR_sigma(isnan(CTD1(i).zHat))=[];
    CTD1(i).zHat(isnan(CTD1(i).zHat))=[];
    
    % Do interpolation based on switch
    switch lower(G_or_L)
        case 'gaussian'
            if length(tmpTR)>2
                [f,gof(i)]=fit(CTD1(i).zHat,tmpTR,'gauss1');
                CTD1(i).c_ij=feval(f,zHat_grid);
                if max(CTD1(i).c_ij)>1 || CTD1(i).c_ij(end)>0.05
                    CTD1(i).c_ij=NaN(length(zHat_grid),1);
                end
                % Create vertical integral
                CTD1(i).c_ij(isnan(CTD1(i).c_ij))=0;
                CTD1(i).I_j=trapz(zHat_grid,CTD1(i).c_ij);
                CTD1(i).r_ij=[CTD1(i).c_ij./CTD1(i).I_j]';
            else
                CTD1(i).I_j=0;
                CTD1(i).c_ij=zeros(length(zHat_grid),1);
                CTD1(i).r_ij=zeros(1,length(zHat_grid));
            end
            
        case 'linear'
            
            if length(tmpTR)>1
                
                if sum(tmpTR)==max(tmpTR)
                    try
                        [f,gof(i)]=fit([CTD1(i).zHat;-15],[tmpTR;0],'gauss1');
                        CTD1(i).c_ij=feval(f,zHat_grid);
                        if max(CTD1(i).c_ij)>1 || CTD1(i).c_ij(end)>0.05
                            CTD1(i).c_ij=NaN(length(zHat_grid),1);
                        end
                    catch
                        CTD1(i).c_ij=NaN(length(zHat_grid),1);
                    end
                else
                
                % interp onto regularly spaced bins
                CTD1(i).c_ij=exp(interp1(CTD1(i).zHat,log(tmpTR),zHat_grid));
                
                end
                
                % Pad with zeros
                CTD1(i).c_ij(isnan(CTD1(i).c_ij))=0;
                
                % Create vertical integral
                CTD1(i).I_j=trapz(zHat_grid,CTD1(i).c_ij);
                CTD1(i).r_ij=[CTD1(i).c_ij./CTD1(i).I_j]';
                
                if max(CTD1(i).r_ij)==1
                    CTD1(i).r_ij=zeros(1,length(zHat_grid));
                end
                
            else
                CTD1(i).I_j=0;
                CTD1(i).c_ij=zeros(length(zHat_grid),1);
                CTD1(i).r_ij=zeros(1,length(zHat_grid));
            end
    end
end


for i=1:length(CTD2)
    % Need unique data for intepolation... add tiny variability
    tmpTR=CTD2(i).TR_SF5;
    tmpSigma=CTD2(i).TR_sigma+(rand(size(CTD2(i).TR_z))/100);
    
    % interpolate these onto CS mean sigma/zHat_grid relationship
    CTD2(i).zStar=interp1(sigCo1,zStar_grid,tmpSigma);
    
    % Change coordinates to distance from injection sigma
    CTD2(i).zHat=CTD2(i).zStar-zStar_grid(sigIdx2);
    CTD2(i).sigHat=CTD2(i).TR_sigma-27.26;
    
    % Remove any data outside of coordinate grid
    tmpTR(isnan(CTD2(i).zStar))=[];
    CTD2(i).zStar(isnan(CTD2(i).zStar))=[];
    CTD2(i).TR_SF5(isnan(CTD2(i).zHat))=[];
    CTD2(i).TR_z(isnan(CTD2(i).zHat))=[];
    CTD2(i).TR_sigma(isnan(CTD2(i).zHat))=[];
    CTD2(i).zHat(isnan(CTD2(i).zHat))=[];
    CTD2(i).sigHat(isnan(CTD2(i).sigHat))=[];
    
    % Do interpolation based on switch
    switch lower(G_or_L)
        case 'gaussian'
            if length(tmpTR)>2
                [f,gof(i)]=fit(CTD2(i).zHat,tmpTR,'gauss1');
                CTD2(i).c_ij=feval(f,zHat_grid);
                if max(CTD2(i).c_ij)>1 || CTD2(i).c_ij(end)>0.05
                    CTD2(i).c_ij=NaN(length(zHat_grid),1);
                end
                % Create vertical integral
                CTD2(i).I_j=trapz(zHat_grid,CTD2(i).c_ij);
                CTD2(i).r_ij=[CTD2(i).c_ij./CTD2(i).I_j]';
                
                % Create mirrored dataset for second cruise- any datapoints within
                % 75m of the injection isopycnal are mirrored
                msk=CTD2(i).zStar>=190 & CTD2(i).zStar<=290;
                CTD2(i).TR_SF5Mir=[tmpTR;tmpTR(msk)];
                CTD2(i).zStarMir=[CTD2(i).zStar;...
                    abs(CTD2(i).zStar(msk)-290)+290];
                CTD2(i).sigMir=[CTD2(i).TR_sigma;...
                    abs(CTD2(i).TR_sigma(msk)-27.26)+27.26];
                CTD2(i).isMir=logical([zeros(size(CTD2(i).zStar));
                    ones(sum(msk),1)]);
                CTD2(i).zHatMir=[abs(CTD2(i).zHat)*-1;...
                    abs(CTD2(i).zHat(msk))];
                
                % interp onto regularly spaced bins
                CTD2(i).c_ijMir=exp(interp1(CTD2(i).zHatMir,...
                    log(CTD2(i).TR_SF5Mir),zHat_grid));
                
                % Zero any NaNs
                CTD2(i).c_ijMir(isnan(CTD2(i).c_ijMir))=0;
                CTD2(i).I_jMir=trapz(zHat_grid(~isnan(CTD2(i).c_ijMir)),...
                    CTD2(i).c_ijMir(~isnan(CTD2(i).c_ijMir)));
                
                CTD2(i).r_ijMir=[CTD2(i).c_ijMir./CTD2(i).I_jMir]';
                
            else
                CTD2(i).I_j=0;
                CTD2(i).I_jMir=0;
                
                CTD2(i).c_ij=zeros(length(zHat_grid),1);
                CTD2(i).c_ijMir=zeros(length(zHat_grid),1);
                
                CTD2(i).r_ij=zeros(1,length(zHat_grid));
                CTD2(i).r_ijMir=zeros(1,length(zHat_grid));
            end
            
        case 'linear'
            if length(tmpTR)>1
                
                if sum(tmpTR)==max(tmpTR)
                    [f,gof(i)]=fit([CTD2(i).zHat;-15],[tmpTR;0],'gauss1');
                    CTD2(i).c_ij=feval(f,zHat_grid);
                    if max(CTD2(i).c_ij)>1 || CTD2(i).c_ij(end)>0.05
                        CTD2(i).c_ij=NaN(length(zHat_grid),1);
                    end
                else
                    % interp onto regularly spaced bins
                    CTD2(i).c_ij=exp(interp1(CTD2(i).zHat,log(tmpTR),zHat_grid));
                end
                
                % Pad with zeros
                CTD2(i).c_ij(isnan(CTD2(i).c_ij))=0;
                
                % Create vertical integral
                CTD2(i).I_j=trapz(zHat_grid,CTD2(i).c_ij);
                CTD2(i).r_ij=[CTD2(i).c_ij./CTD2(i).I_j]';
                
                % Create mirrored dataset for second cruise- any datapoints within
                % 75m of the injection isopycnal are mirrored
                msk=CTD2(i).zStar>=190 & CTD2(i).zStar<=290;
                CTD2(i).TR_SF5Mir=[tmpTR;tmpTR(msk)];
                CTD2(i).zStarMir=[CTD2(i).zStar;...
                    abs(CTD2(i).zStar(msk)-290)+290];
                CTD2(i).sigMir=[CTD2(i).TR_sigma;...
                    abs(CTD2(i).TR_sigma(msk)-27.26)+27.26];
                CTD2(i).isMir=logical([zeros(size(CTD2(i).zStar));
                    ones(sum(msk),1)]);
                CTD2(i).zHatMir=[abs(CTD2(i).zHat)*-1;...
                    abs(CTD2(i).zHat(msk))];
                
                % interp onto regularly spaced bins
                CTD2(i).c_ijMir=exp(interp1(CTD2(i).zHatMir,...
                    log(CTD2(i).TR_SF5Mir),zHat_grid));
                
                % Zero any NaNs
                CTD2(i).c_ijMir(isnan(CTD2(i).c_ijMir))=0;
                CTD2(i).I_jMir=trapz(zHat_grid(~isnan(CTD2(i).c_ijMir)),...
                    CTD2(i).c_ijMir(~isnan(CTD2(i).c_ijMir)));
                
                CTD2(i).r_ijMir=[CTD2(i).c_ijMir./CTD2(i).I_jMir]';
                
            else
                CTD2(i).I_j=0;
                CTD2(i).I_jMir=0;
                
                CTD2(i).c_ij=zeros(length(zHat_grid),1);
                CTD2(i).c_ijMir=zeros(length(zHat_grid),1);
                
                CTD2(i).r_ij=zeros(1,length(zHat_grid));
                CTD2(i).r_ijMir=zeros(1,length(zHat_grid));
            end
    end
end


%% Ledwell and Watson equations 3 and 4 weighted averages for cruise
c_overbar1=mean(horzcat(CTD1.c_ij),2,'omitnan');
for i=1:length(CTD1)
    CTD1(i).w_j=CTD1(i).I_j/sum(vertcat(CTD1.I_j));
    for ii=1:length(CTD1(i).TR_SF5)
        [~,idx]=min(abs(CTD1(i).zHat(ii)-zHat_grid));
        CTD1(i).n_ij(ii)=CTD1(i).TR_SF5(ii)./c_overbar1(idx);
    end
end

c_overbar2=mean(horzcat(CTD2.c_ij),2,'omitnan');
for i=1:length(CTD2)
    CTD2(i).w_j=CTD2(i).I_j/sum(vertcat(CTD2.I_j));
    for ii=1:length(CTD2(i).TR_SF5)
        [~,idx]=min(abs(CTD2(i).zHat(ii)-zHat_grid));
        CTD2(i).n_ij(ii)=CTD2(i).TR_SF5(ii)./c_overbar2(idx);
    end
end

for i=1:length(CTD2)
    CTD2(i).w_jMir=CTD2(i).I_jMir/sum(vertcat(CTD2.I_jMir));
end

allwj1=vertcat(CTD1.w_j);
allwj2=vertcat(CTD2.w_j);
allwj2Mir=vertcat(CTD2.w_jMir);

allrij1=vertcat(CTD1.r_ij);
allrij2=vertcat(CTD2.r_ij);
allrij2Mir=vertcat(CTD2.r_ijMir);

ri1=sum(allwj1.*allrij1,'omitnan');
ri2=sum(allwj2.*allrij2,'omitnan');
ri2Mir=sum(allwj2Mir.*allrij2Mir,'omitnan');

%% Bootstrap some errors on <r> profile
I1=bootrnd(size(allrij1,1),1001);
I2=bootrnd(size(allrij2,1),1001);

errAll1=NaN(length(zHat_grid),1000);
errAll2=NaN(length(zHat_grid),1000);
errAll2Mir=NaN(length(zHat_grid),1000);

for i=2:length(I1)
    errAll1(:,i-1)=std(allrij1(I1(1:ceil(size(allrij1,1)/2),i),:)...
        .*allwj1(I1(1:ceil(size(allrij1,1)/2),i)),0,'omitnan');
    
    errAll2(:,i-1)=std(allrij2(I2(1:ceil(size(allrij2,1)/2),i),:)...
        .*allwj2(I2(1:ceil(size(allrij2,1)/2),i)),0,'omitnan');
    
    errAll2Mir(:,i-1)=std(allrij2Mir(I2(1:ceil(size(allrij2Mir,1)/2),i),:)...
        .*allwj2Mir(I2(1:ceil(size(allrij2Mir,1)/2),i)),0,'omitnan');
end

r_errLim1=mean(errAll1,2);
r_errLim2=mean(errAll2,2);
r_errLim2Mir=mean(errAll2Mir,2);

%% Create interpolated Ij fields

% Load bathy and pull out station locs
bathy=load('GoSL_Zhires.mat');
bathy.Z(bathy.Z==-9999)=NaN;
% Decimate bathy to speed things up
dec=10;
bathy.lon=bathy.lon(1:dec:end,1:dec:end);
bathy.lat=bathy.lat(1:dec:end,1:dec:end);
bathy.Z=bathy.Z(1:dec:end,1:dec:end);

% create grids with 200 m (?) spacing
[xx,yy]=meshgrid(-69.5:0.02:-57,46:0.02:51.5);
lat_lim=[46.75 51.5];
lon_lim=[-69.5 -57.5];
m_proj('UTM','lon',lon_lim,'lat',lat_lim);
[xm,ym]=m_ll2xy(xx,yy);
[x,y]=m_ll2xy(stlon1,stlat1);

% Create nan mask
msk1=~isnan(vertcat(CTD1.I_j));
z=vertcat(CTD1.I_j);

% Try kriging
dday = variogram([x(msk1) y(msk1)],z(msk1),'plotit',false);
[a,~,n,vstruct] = variogramfit(dday.distance,dday.val,200e3,4e-8,dday.num,'plotit',false);
[zi,zivar] = kriging(vstruct,x(msk1),y(msk1),z(msk1),xm,ym);

[elev,Zlon,Zlat]=m_etopo2([lon_lim lat_lim]);
elevI=interp2(Zlon,Zlat,elev,xx,yy);
zi(elevI>-80)=NaN;

% zi(xx<62 & yy<47.8)=NaN;
zi(zi<0)=0;

%%%%%%
[x,y]=m_ll2xy(stlon2,stlat2);

% Create nan mask
msk2=~isnan(vertcat(CTD2.I_j));
z=vertcat(CTD2.I_j);

% Try kriging
dday = variogram([x(msk2) y(msk2)],z(msk2),'plotit',false);
[a,~,n,vstruct] = variogramfit(dday.distance,dday.val,200e3,4e-8,dday.num,'plotit',false);
[zi2,zivar] = kriging(vstruct,x(msk2),y(msk2),z(msk2),xm,ym);

[elev,Zlon,Zlat]=m_etopo2([lon_lim lat_lim]);
elevI=interp2(Zlon,Zlat,elev,xx,yy);
zi2(elevI>-80)=NaN;

% zi2(xx<62 & yy<47.8)=NaN;
zi2(zi2<0)=0;

% Remove data outside measurement region
as=alphaShape([stlon1;stlon2(stlat2>47.75 & stlat2<max(stlat1) & stlon2>-67.75)],...
    [stlat1;stlat2(stlat2>47.75 & stlat2<max(stlat1) & stlon2>-67.75)]...
    ,1,'HoleThreshold',15);
zi(~inShape(as,xx,yy))=NaN;
as=alphaShape([stlon1;stlon2],[stlat1;stlat2],1,'HoleThreshold',15);
zi2(~inShape(as,xx,yy))=NaN;

%%
save('/Users/samst/Library/CloudStorage/Dropbox/UBC/TReX_paper/code/letter/data20230707.mat');

