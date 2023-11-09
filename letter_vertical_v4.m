%% Vertical mixing analysis
addpath(genpath('/Users/samst/Dropbox/UBC/TReX_deepII/'));
addpath(genpath('/Users/samst/Dropbox/UBC/TReX_paper/'));
addpath(genpath('/Users/samst/Dropbox/UBC/GSW/'));
addpath(genpath('/Users/samst/Dropbox/UBC/m_map/'));
addpath(genpath('/Users/samst/Dropbox/UBC/Misc/'));

clear

load data20230707.mat
model10=load('model10_z.mat');

% Create model fields
model10.Cp=squeeze(sum(sum(model10.CDzi,'omitnan')));

%% Bootstrap some errors on <r> profile
I1=bootrnd(size(allrij1,1),1001);
I2=bootrnd(size(allrij2,1),1001);

errAll1=NaN(length(zHat_grid),1000); errAll2=NaN(length(zHat_grid),1000);
errAll2Mir=NaN(length(zHat_grid),1000);

for i=2:length(I1)
    errAll1(:,i-1)=std(allrij1(I1(1:ceil(size(allrij1,1)/2),i),:)...
        .*allwj1(I1(1:ceil(size(allrij1,1)/2),i)),0,'omitnan');
    
    errAll2(:,i-1)=std(allrij2(I2(1:ceil(size(allrij2,1)/2),i),:)...
        .*allwj2(I2(1:ceil(size(allrij2,1)/2),i)),0,'omitnan');
    
    errAll2Mir(:,i-1)=std(allrij2Mir(I2(1:ceil(size(allrij2Mir,1)/2),i),:)...
        .*allwj2Mir(I2(1:ceil(size(allrij2Mir,1)/2),i)),0,'omitnan');
end

errLim1=mean(errAll1,2).*2;
errLim2=mean(errAll2,2).*2;
errLim2Mir=mean(errAll2Mir,2).*2;

%% Calculate w

t=[(CTD1(1).mtime(1)-datenum(2021,10,24)) ...
    (CTD2(1).mtime(1)-datenum(2021,10,24))];  % day

[~,idx1]=max(ri1);[~,idx2]=max(ri2);

wa=[zHat_grid(idx1) zHat_grid(idx2)]./(t*86400); % m/s

% w is constant between cruises, deepening at ~17 m/yr
w=mean(wa);

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

%% Advection-diffusion model
M0=[71 60]; % g

A0=[Az_star1(1) Az_star2(1)];% m^2

%%%%%
L=200; % m
kappa=[0.05e-5:0.1e-5:10e-5].*86400; % m^2/day
U=-25/365/86400:0.5/365/86400:25/365/86400; % m/s

ri_model1=model10.Cp/sum(model10.Cp);

% This weighting emphasises shallow regions of the distribution and removes
% the influence of depths with no measurements
h1=histcounts(vertcat(CTD1.zHat),zHat_grid);
h2=histcounts(vertcat(CTD2.zHat),zHat_grid);
wey=[ones(size(h1')) ones(size(h2'))];
wey(h1==0,1)=0;wey(h2==0,2)=0;
wey=[wey;[0 0]];
wey(zHat_grid>0,2)=0;

idxs=NaN(2,length(U)*length(kappa));

% run through different kappas and find the best fit
cost1=NaN(length(kappa)*length(U),1);
cost2=NaN(length(kappa)*length(U),1);
cost3=NaN(length(kappa)*length(U),1);
cost1E=NaN(length(kappa)*length(U),2);
cost2E=NaN(length(kappa)*length(U),2);

count=0;
for i=1:length(kappa)
    for ii=1:length(U)
        count=count+1;
        
        w=U(ii);
        W=(w*86400)+kappa(i)/L; % m/day
        
        K(count).c=M0./A0./sqrt(4*pi*kappa(i)*t).*exp(-(zHat_grid-W*t).^2./(4*kappa(i)*t) + (w*86400)*t/L )...
            ./1e3.*1e15./mm_sf5cf3; % fmol/l
        
        % convert to fmol/kg using average CS density profile (very similar in
        % deeper waters)
        K(count).c=K(count).c.*((26+1000)/1000); %fmol/kg
        
        % Find Least square difference between ris and model
        Ij=trapz(zHat_grid,K(count).c(:,1));
        K(count).rM1=K(count).c(:,1)./Ij;%.*trapz(z,ri1);
        
        % Cruise 1 method
        cost1(count)=trapz(zHat_grid,wey(:,1).*(K(count).rM1-...
            ri1').^2);
        cost1E(count,1)=trapz(zHat_grid,wey(:,1).*(K(count).rM1-...
            ri1'-errLim1).^2);
        cost1E(count,2)=trapz(zHat_grid,wey(:,1).*(K(count).rM1-...
            ri1'+errLim1).^2);
        
        % For PT model
        costM(count)=trapz(zHat_grid,(K(count).rM1-ri_model1).^2);
        
        %Cruise 2
        Ij=trapz(zHat_grid,K(count).c(:,2));
        K(count).rM2=K(count).c(:,2)./Ij;%.*trapz(z,ri2);
        
        % Normal
        cost2(count)=trapz(zHat_grid,wey(:,2).*(K(count).rM2-...
            ri2').^2);
         cost2E(count,1)=trapz(zHat_grid,wey(:,1).*(K(count).rM2-...
            ri2'-errLim2).^2);
        cost2E(count,2)=trapz(zHat_grid,wey(:,1).*(K(count).rM2-...
            ri2'+errLim2).^2);
        
        % Mirrored fitting
        cost3(count)=trapz(zHat_grid,(K(count).rM2-ri2Mir').^2);
        
        idxs(:,count)=[i;ii];
        
    end
end

% find minimum cost
[~,Cidx1]=min(cost1);
[~,Cidx1_1]=min(cost1E(:,1));
[~,Cidx1_2]=min(cost1E(:,2));
[~,Cidx2]=min(cost2);
[~,Cidx2_1]=min(cost2E(:,1));
[~,Cidx2_2]=min(cost2E(:,2));
[~,Cidx3]=min(cost3);
[~,CidxM]=min(costM);

clc
fprintf('Leg 1 1D model: k=%2.1e (%2.1e-%2.1e) m^2 s^{-1}, w=%2.1f (%2.1f-%2.1f) m yr^{-1}\n',...
    kappa(idxs(1,Cidx1))/86400,kappa(idxs(1,Cidx1_1))/86400,kappa(idxs(1,Cidx1_2))/86400,...
    U(idxs(2,Cidx1))*365*86400,U(idxs(2,Cidx1_1))*365*86400,U(idxs(2,Cidx1_2))*365*86400);
fprintf('Leg 1 1D model (CIOPS-E): k=%2.1e m^2 s^{-1}, w=%2.1f m yr^{-1}\n',...
    kappa(idxs(1,CidxM))/86400,U(idxs(2,CidxM))*365*86400);
fprintf('Leg 2 1D model: k=%2.1e (%2.1e-%2.1e) m^2 s^{-1}, w=%2.1f (%2.1f-%2.1f) m yr^{-1}\n',...
    kappa(idxs(1,Cidx2))/86400,kappa(idxs(1,Cidx2_1))/86400,kappa(idxs(1,Cidx2_2))/86400,...
    U(idxs(2,Cidx2))*365*86400,U(idxs(2,Cidx2_1))*365*86400,U(idxs(2,Cidx2_2))*365*86400);

% Calculate the diapycnal diffusivity
d_rho_d_zStar = gradient(sigCo1, zStar_grid);
kappa_rho=(kappa(idxs(1,Cidx1))/86400)*mean(d_rho_d_zStar(200:375).^2);

kappa_L1=(kappa(idxs(1,Cidx1))*365)/L
kappa_L2=(kappa(idxs(1,Cidx2))*365)/L


%% Bootstrap errors
kerr1=NaN(size(I1,2),1); Uerr1=NaN(size(I1,2),1);
kerr2=NaN(size(I2,2),1); Uerr2=NaN(size(I2,2),1);

hw=waitbar(0,'Calculating...');
tic

allIj1=vertcat(CTD1.I_j);
allIj2=vertcat(CTD2.I_j);

for j=2:size(I1,2)
    tmpWj1=allIj1(I1(:,j))/sum(allIj1(I1(:,j)));
    r1_tmp=sum(allrij1(I1(:,j),:).*tmpWj1,1,'omitnan');
    
    tmpWj2=allIj2(I2(:,j))/sum(allIj2(I2(:,j)));
    r2_tmp=sum(allrij2(I2(:,j),:).*tmpWj2,1,'omitnan');
    
    idxsE=NaN(2,length(U)*length(kappa));
    count=0;
    cost1E=NaN(1,length(U)*length(kappa));
    cost2E=NaN(1,length(U)*length(kappa));
    
    for i=1:length(K)
        count=count+1;
        % Cruise 1
        cost1E(count)=trapz(zHat_grid,wey(:,1).*(K(i).rM1-...
            r1_tmp').^2);
        
        %Cruise 2
        cost2E(count)=trapz(zHat_grid,wey(:,2).*(K(i).rM2-...
            r2_tmp').^2);
    end
    [~,idxC1E]=min(cost1E); [~,idxC2E]=min(cost2E);
    
%     figure;hold on;
%     xlim([0 0.03])
%     plot(r1_tmp,zHat_grid,'r--');
%     plot(K(idxC1E).rM1,zHat_grid,'r');
%     plot(r2_tmp,zHat_grid,'b--');
%     plot(K(idxC2E).rM2,zHat_grid,'b');
%     pause
%     
%     figure;plot(cost1E)
%     hold on
%     scatter(idxC1E,min(cost1E),20,'filled');
%     pause
    
    kerr1(j,1)=kappa(idxs(1,idxC1E));
    Uerr1(j,1)=U(idxs(2,idxC1E));
    kerr2(j,1)=kappa(idxs(1,idxC2E));
    Uerr2(j,1)=U(idxs(2,idxC2E));
    
    waitbar(j/size(I1,2),hw);
    
end
toc
delete(hw)

save('/Users/samst/Library/CloudStorage/Dropbox/UBC/TReX_paper/data/1DbootZ.mat',...
    'kerr1', 'Uerr1', 'kerr2', 'Uerr2');

load 1DbootZ.mat

% From 1000 subsets
fprintf('\nLeg 1 (2 STD): k=%2.1e m^2 s^{-1}, U=%4.3f m yr^{-1}\n',...
    nanstd(kerr1)*2/86400,nanstd(Uerr1)*2*365*86400);
fprintf('\nLeg 2 (2 STD): k=%2.1e m^2 s^{-1}, U=%4.3f m yr^{-1}\n',...
    nanstd(kerr2)/86400*2,nanstd(Uerr2)*2*365*86400);

%% What about deep vs shallow? Cruise 1
deepIj1=[]; shallowIj1=[];
deepRij1=[]; shallowRij1=[];

for i=1:length(CTD1)
    if stBath1(i)<-300 && CTD1(i).I_j>0
        deepIj1=[deepIj1;CTD1(i).I_j];
        deepRij1=[deepRij1;CTD1(i).r_ij];
        
    elseif stBath1(i)>-300 && CTD1(i).I_j>0
        shallowIj1=[shallowIj1;CTD1(i).I_j];
        shallowRij1=[shallowRij1;CTD1(i).r_ij];
        
    end
end

% Do normalization for different depths
shallow_wj1=shallowIj1./...
    sum(shallowIj1);
shallow_ri1=sum(shallow_wj1.*shallowRij1,'omitnan');

deep_wj1=deepIj1./...
    sum(deepIj1);
deep_ri1=sum(deep_wj1.*deepRij1,'omitnan');

% Bootstrap some errors
I1=bootrnd(size(shallowRij1,1),1001);
I2=bootrnd(size(deepRij1,1),1001);

errShallow1=NaN(length(zHat_grid),1000);
errDeep1=NaN(length(zHat_grid),1000);

for i=2:length(I1)
    errShallow1(:,i-1)=std(shallowRij1(I1(1:ceil(size(shallowRij1,1)/2),i),:)...
        .*shallow_wj1(I1(1:ceil(size(shallowRij1,1)/2),i)),0,'omitnan');
    
    errDeep1(:,i-1)=std(deepRij1(I2(1:ceil(size(deepRij1,1)/2),i),:)...
        .*deep_wj1(I2(1:ceil(size(deepRij1,1)/2),i)),0,'omitnan');
end

errShallow1=mean(errShallow1,2).*2;
errDeep1=mean(errDeep1,2).*2;

% cost fitting
cost1=trapz(zHat_grid,([K.rM1]-shallow_ri1').^2);
cost2=trapz(zHat_grid,([K.rM1]-deep_ri1').^2);
[~,CidxDeep1]=min(cost2);
[~,CidxShallow1]=min(cost1);

fprintf('Leg 1 shallow:\n k=%2.1e m^2 s^{-1}, w=%2.1f m yr^{-1}\n',...
    kappa(idxs(1,CidxShallow1))/86400,U(idxs(2,CidxShallow1))*365*86400);
fprintf('Leg 1 deep:\n k=%2.1e m^2 s^{-1}, w=%2.1f m yr^{-1}\n',...
    kappa(idxs(1,CidxDeep1))/86400,U(idxs(2,CidxDeep1))*365*86400)

kerr1=NaN(size(I1,2),1); Uerr1=NaN(size(I1,2),1);
kerr2=NaN(size(I2,2),1); Uerr2=NaN(size(I2,2),1);

hw=waitbar(0,'Calculating...');
tic

for j=2:size(I1,2)
    tmpWjS=shallowIj1(I1(:,j))/sum(shallowIj1(I1(:,j)));
    rS_tmp=sum(shallowRij1(I1(:,j),:).*tmpWjS,1,'omitnan');
    
    tmpWjD=deepIj1(I2(:,j))/sum(deepIj1(I2(:,j)));
    rD_tmp=sum(deepRij1(I2(:,j),:).*tmpWjD,1,'omitnan');
    
    idxsE=NaN(2,length(U)*length(kappa));
    count=0;
    cost1E=NaN(1,length(U)*length(kappa));
    cost2E=NaN(1,length(U)*length(kappa));
    
    for i=1:length(K)
        count=count+1;
        % Cruise 1
        cost1E(count)=trapz(zHat_grid,wey(:,1).*(K(i).rM1-...
            rS_tmp').^2);
        
        %Cruise 2
        cost2E(count)=trapz(zHat_grid,wey(:,2).*(K(i).rM2-...
            rD_tmp').^2);
    end
    [~,idxC1E]=min(cost1E); [~,idxC2E]=min(cost2E);
    
%     figure;hold on;
%     xlim([0 0.03])
%     plot(r1_tmp,zHat_grid,'r--');
%     plot(K(idxC1E).rM1,zHat_grid,'r');
%     plot(r2_tmp,zHat_grid,'b--');
%     plot(K(idxC2E).rM2,zHat_grid,'b');
%     pause
%     
%     figure;plot(cost1E)
%     hold on
%     scatter(idxC1E,min(cost1E),20,'filled');
%     pause
%     
    kerr1(j,1)=kappa(idxs(1,idxC1E));
    Uerr1(j,1)=U(idxs(2,idxC1E));
    kerr2(j,1)=kappa(idxs(1,idxC2E));
    Uerr2(j,1)=U(idxs(2,idxC2E));
    
    waitbar(j/size(I1,2),hw);
    
end
toc
delete(hw)

% From 1000 subsets
fprintf('\nLeg 1 boundary error (2 STD): k=%2.1e m^2 s^{-1}, U=%4.3f m yr^{-1}\n',...
    nanstd(kerr1)*2/86400,nanstd(Uerr1)*2*365*86400);
fprintf('\nLeg 1 interior error (2 STD): k=%2.1e m^2 s^{-1}, U=%4.3f m yr^{-1}\n',...
    nanstd(kerr2)*2/86400,nanstd(Uerr2)*2*365*86400);

%%

% CIOPS-E
mTracer=load('conc_20220604_K10');
bathyC=griddata(bathy.lon,bathy.lat,bathy.Z,mTracer.lon,mTracer.lat);
msk=find(bathyC<300);

% Do normalization for different depths
shallow_wj1=sum(model10.CDzi(msk),3,'omitnan')./...
    sum(model10.CDzi(msk),1:3,'omitnan');

modelS=NaN(401,length(msk));
shallowRij=NaN(401,length(msk));
for i=1:length(msk)
    [row,col] = ind2sub(size(model10.CDzi),msk(i));
    modelS(:,i)=model10.CDzi(row,col,:);
    shallowRij(:,i)=modelS(:,i)./sum(modelS(:,i),'omitnan');
end
shallow_ri1=sum((shallow_wj1'.*shallowRij),2,'omitnan');

msk=find(bathyC>300);
% Do normalization for different depths
deep_wj1=sum(model10.CDzi(msk),3,'omitnan')./...
    sum(model10.CDzi(msk),1:3,'omitnan');
modelD=NaN(401,length(msk));
deepRij=NaN(401,length(msk));
for i=1:length(msk)
    [row,col] = ind2sub(size(model10.CDzi),msk(i));
    modelD(:,i)=model10.CDzi(row,col,:);
    deepRij(:,i)=modelD(:,i)./sum(modelD(:,i),'omitnan');
end
deep_ri1=sum((deep_wj1'.*deepRij),2,'omitnan');
%%
% Leg 2
deepIj2=[]; shallowIj2=[];
deepRij2=[]; shallowRij2=[];

for i=1:length(CTD2)
    if stBath2(i)<-300 && CTD2(i).I_j>0
        deepIj2=[deepIj2;CTD2(i).I_j];
        deepRij2=[deepRij2;CTD2(i).r_ij];
        
    elseif stBath2(i)>-300 && CTD2(i).I_j>0
        shallowIj2=[shallowIj2;CTD2(i).I_j];
        shallowRij2=[shallowRij2;CTD2(i).r_ij];
        
    end
end

% Do normalization for different depths
shallow_wj2=shallowIj2./...
    sum(shallowIj2);
shallow_ri2=sum(shallow_wj2.*shallowRij2,'omitnan');

deep_wj2=deepIj2./...
    sum(deepIj2);
deep_ri2=sum(deep_wj2.*deepRij2,'omitnan');

% Bootstrap some errors
I1=bootrnd(size(shallowRij2,1),1001);
I2=bootrnd(size(deepRij2,1),1001);

errShallow2=NaN(length(zHat_grid),1000);
errDeep2=NaN(length(zHat_grid),1000);

for i=2:length(I1)
    errShallow2(:,i-1)=std(shallowRij2(I1(1:ceil(size(shallowRij2,1)/2),i),:)...
        .*shallow_wj2(I1(1:ceil(size(shallowRij2,1)/2),i)),0,'omitnan');
    
    errDeep2(:,i-1)=std(deepRij2(I2(1:ceil(size(deepRij2,1)/2),i),:)...
        .*deep_wj2(I2(1:ceil(size(deepRij2,1)/2),i)),0,'omitnan');
end

errShallow2=mean(errShallow2,2).*2;
errDeep2=mean(errDeep2,2).*2;

% cost fitting
cost1=trapz(zHat_grid,wey(:,2).*([K.rM2]-shallow_ri2').^2);
[~,CidxShallow2]=min(cost1);

cost2=trapz(zHat_grid,wey(:,2).*([K.rM2]-deep_ri2').^2);
[~,CidxDeep2]=min(cost2);

fprintf('\nLeg 2 shallow:\n k=%2.1e m^2 s^{-1}, w=%2.1f m yr^{-1}\n',...
    kappa(idxs(1,CidxShallow2))/86400,U(idxs(2,CidxShallow2))*365*86400);
fprintf('Leg 2 deep:\n k=%2.1e m^2 s^{-1}, w=%2.1f m yr^{-1}\n',...
    kappa(idxs(1,CidxDeep2))/86400,U(idxs(2,CidxDeep2))*365*86400)

fprintf('\n\nmean shallow:\n k=%2.1e m^2 s^{-1}, w=%2.1f m yr^{-1}\n',...
    mean([kappa(idxs(1,CidxShallow2)) kappa(idxs(1,CidxShallow1))])/86400,...
    mean([U(idxs(2,CidxShallow2)) U(idxs(2,CidxShallow1))])*365*86400);

fprintf('\nmean deep:\n k=%2.1e m^2 s^{-1}, w=%2.1f m yr^{-1}\n',...
    mean([kappa(idxs(1,CidxDeep2)) kappa(idxs(1,CidxDeep1))])/86400,...
    mean([U(idxs(2,CidxDeep2)) U(idxs(2,CidxDeep1))])*365*86400);

fprintf('\n Deep subset is %2.1f pcnt of profiles (both legs)\n',...
    (length([deepIj1;deepIj2])/...
    length([deepIj1;deepIj2;shallowIj1;shallowIj2]))*100);


kerr1=NaN(size(I1,2),1); Uerr1=NaN(size(I1,2),1);
kerr2=NaN(size(I2,2),1); Uerr2=NaN(size(I2,2),1);

hw=waitbar(0,'Calculating...');
tic

for j=2:size(I1,2)
    tmpWjS=shallowIj2(I1(:,j))/sum(shallowIj2(I1(:,j)));
    rS_tmp=sum(shallowRij2(I1(:,j),:).*tmpWjS,1,'omitnan');
    
    tmpWjD=deepIj2(I2(:,j))/sum(deepIj2(I2(:,j)));
    rD_tmp=sum(deepRij2(I2(:,j),:).*tmpWjD,1,'omitnan');
    
    idxsE=NaN(2,length(U)*length(kappa));
    count=0;
    cost1E=NaN(1,length(U)*length(kappa));
    cost2E=NaN(1,length(U)*length(kappa));
    
    for i=1:length(K)
        count=count+1;
        % Cruise 1
        cost1E(count)=trapz(zHat_grid,wey(:,1).*(K(i).rM1-...
            rS_tmp').^2);
        
        %Cruise 2
        cost2E(count)=trapz(zHat_grid,wey(:,2).*(K(i).rM2-...
            rD_tmp').^2);
    end
    [~,idxC1E]=min(cost1E); [~,idxC2E]=min(cost2E);
    
%     figure;hold on;
%     xlim([0 0.03])
%     plot(r1_tmp,zHat_grid,'r--');
%     plot(K(idxC1E).rM1,zHat_grid,'r');
%     plot(r2_tmp,zHat_grid,'b--');
%     plot(K(idxC2E).rM2,zHat_grid,'b');
%     pause
%     
%     figure;plot(cost1E)
%     hold on
%     scatter(idxC1E,min(cost1E),20,'filled');
%     pause
%     
    kerr1(j,1)=kappa(idxs(1,idxC1E));
    Uerr1(j,1)=U(idxs(2,idxC1E));
    kerr2(j,1)=kappa(idxs(1,idxC2E));
    Uerr2(j,1)=U(idxs(2,idxC2E));
    
    waitbar(j/size(I1,2),hw);
    
end
toc
delete(hw)

% From 1000 subsets
fprintf('\nLeg 2 boundary error (2 STD): k=%2.1e m^2 s^{-1}, U=%4.3f m yr^{-1}\n',...
    nanstd(kerr1)*2/86400,nanstd(Uerr1)*2*365*86400);
fprintf('\nLeg 2 interior error (2 STD): k=%2.1e m^2 s^{-1}, U=%4.3f m yr^{-1}\n',...
    nanstd(kerr2)*2/86400,nanstd(Uerr2)*2*365*86400);

%% Regional fitting
f1=figure('units','centimeters','outerposition',...
    [0.01 0.01 19 23],'color','w');
col=[255 214 140]/255; % YELLOW!

ax1=axes('position',[0.125 0.3 .75 .75]); hold on;
lat_lim=[46 51];
lon_lim=[-67.5 -57.5];
m_proj('equidistant','lon',lon_lim,'lat',lat_lim); hold on;
[c,h]=m_contour(bathy.lon,bathy.lat,bathy.Z,0:-100:-500,'color',...
    [0.95 0.95 0.95],'linewidth',0.5);
clabel(c,h,'LabelSpacing',500,'color',[0.95 0.95 0.95],'fontsize',6);
m_gshhs_f('patch',col,'edgecolor',[0.3 0.3 0.3],'linewidth',0.2); % CHANGE
m_grid('linestyle','none','linewidth',0.1,'tickdir','out',...
    'xaxisloc','top','yaxisloc','left','fontsize',6);
m_text(-67,50.75,'a)','fontsize',8);

%%%%%%%%%%%%%%%%%%%%%%% Define shapes
% Define WLA shape
WL.lon=[-66.4852,-65.8573,-65.1713,-64.4853,-63.9970,-63.4506,-64.2761,...
    -65.4969,-66.1712,-66.9735];
WL.lat=[49.3725,49.2800,49.3032,49.1568,48.8794,49.3802,49.7808,49.8964,...
    50.0813,49.7038];
WL.station_bound=alphaShape(WL.lon',WL.lat',1,'HoleThreshold',15);
[l,ll]=m_ll2xy(WL.lon,WL.lat);
tmpshp=alphaShape(l',ll',1,'HoleThreshold',15);
WL.AS=plot(tmpshp,'facealpha',0.1,'edgecolor','none','facecolor','g');
% find the best fit
WL.idx1=inShape(WL.station_bound,stlon1,stlat1);
WL.z1=vertcat(CTD1(WL.idx1).zStar);WL.tr1=vertcat(CTD1(WL.idx1).TR_SF5);
WL.idx2=inShape(WL.station_bound,stlon2,stlat2);
WL.z2=vertcat(CTD2(WL.idx2).zStar);WL.tr2=vertcat(CTD2(WL.idx2).TR_SF5);
WL.cost1=trapz(zHat_grid,([K.c]-nanmean([CTD1(WL.idx1).c_ij],2)).^2);
WL.cost2=trapz(zHat_grid,([K.c]-nanmean([CTD2(WL.idx2).c_ij],2)).^2);
[~,WL.Ridx1]=min(cost1); [~,WL.Ridx2]=min(cost2);

% Define ELA shape
EL.lon=[-63.4622,-63.9738,-63.8110,-62.8460,-62.5669,-63.4506];
EL.lat=[49.3802,48.8794,48.5636,48.4711,49.1183,49.3725];
EL.station_bound=alphaShape(EL.lon',EL.lat',1,'HoleThreshold',15);
[l,ll]=m_ll2xy(EL.lon,EL.lat);
tmpshp=alphaShape(l',ll',1,'HoleThreshold',15);
EL.AS=plot(tmpshp,'facealpha',0.1,'edgecolor','none','facecolor',rgb_x('dark teal'));
% find the best fit
EL.idx1=inShape(EL.station_bound,stlon1,stlat1);
EL.z1=vertcat(CTD1(EL.idx1).zStar);EL.tr1=vertcat(CTD1(EL.idx1).TR_SF5);
EL.idx2=inShape(EL.station_bound,stlon2,stlat2);
EL.z2=vertcat(CTD2(EL.idx2).zStar);EL.tr2=vertcat(CTD2(EL.idx2).TR_SF5);
EL.cost1=trapz(zHat_grid,([K.c]-nanmean([CTD1(EL.idx1).c_ij],2)).^2);
EL.cost2=trapz(zHat_grid,([K.c]-nanmean([CTD2(EL.idx2).c_ij],2)).^2);
[~,EL.Ridx1]=min(cost1); [~,EL.Ridx2]=min(cost2);

% Define OL shape
OL.lon=[-62.2530,-62.5786,-61.9507,-61.4275,-60.6601,...
    -61.4159,-62.5786,-62.8576];
OL.lat=[49.1106+0.02,48.2246+0.1,48.2246,48.0859,48.7870,...
    48.9488,49.1106,48.4634];
OL.station_bound=alphaShape(OL.lon',OL.lat',1,'HoleThreshold',15);
[l,ll]=m_ll2xy(OL.lon,OL.lat);
tmpshp=alphaShape(l',ll',1,'HoleThreshold',15);
OL.AS=plot(tmpshp,'facealpha',0.1,'facecolor','r','edgecolor','none');
% find the best fit
OL.idx1=inShape(OL.station_bound,stlon1,stlat1);
OL.z1=vertcat(CTD1(OL.idx1).zStar);OL.tr1=vertcat(CTD1(OL.idx1).TR_SF5);
OL.idx2=inShape(OL.station_bound,stlon2,stlat2);
OL.z2=vertcat(CTD2(OL.idx2).zStar);OL.tr2=vertcat(CTD2(OL.idx2).TR_SF5);
OL.cost1=trapz(zHat_grid,([K.c]-nanmean([CTD1(OL.idx1).c_ij],2)).^2);
OL.cost2=trapz(zHat_grid,([K.c]-nanmean([CTD2(OL.idx2).c_ij],2)).^2);
[~,OL.Ridx1]=min(cost1); [~,OL.Ridx2]=min(cost2);

% Define CS shape
CS.lon=[-60.6601,-61.4159,-60.7996,-60.1602,-59.3463,-59.3928,-59.6253,-57.5674,...
    -58.7998,-58.0441,-59.8579];
CS.lat=[48.7870,48.0859,47.8163,47.2461,47.5851,48.1629,48.3401,46.5913,...
    46.0289,47.2153,46.6529];
CS.station_bound=alphaShape(CS.lon',CS.lat',1,'HoleThreshold',15);
[l,ll]=m_ll2xy(CS.lon,CS.lat);
tmpshp=alphaShape(l',ll',1,'HoleThreshold',15);
CS.AS=plot(tmpshp,'facealpha',0.15,'facecolor',rgb_x('light gold'),'edgecolor','none');
% find the best fit
CS.idx1=inShape(CS.station_bound,stlon1,stlat1);
CS.z1=vertcat(CTD1(CS.idx1).zStar);CS.tr1=vertcat(CTD1(CS.idx1).TR_SF5);
CS.idx2=inShape(CS.station_bound,stlon2,stlat2);
CS.z2=vertcat(CTD2(CS.idx2).zStar);CS.tr2=vertcat(CTD2(CS.idx2).TR_SF5);
CS.cost1=trapz(zHat_grid,([K.c]-nanmean([CTD1(CS.idx1).c_ij],2)).^2);
CS.cost2=trapz(zHat_grid,([K.c]-nanmean([CTD2(CS.idx2).c_ij],2)).^2);
[~,CS.Ridx1]=min(cost1); [~,CS.Ridx2]=min(cost2);

% Define BL shape
BL.lon=[-58.3580,-59.2998,-59.7416-0.2,-58.7766,-57.5209,-59.8230,-57.5325];
BL.lat=[50.6822,50.1583,49.6653-0.1,49.5343-0.25,50.5281,49.9580,50.8286];
BL.station_bound=alphaShape(BL.lon',BL.lat',1,'HoleThreshold',15);
[l,ll]=m_ll2xy(BL.lon,BL.lat);
tmpshp=alphaShape(l',ll',1,'HoleThreshold',15);
BL.AS=plot(tmpshp,'facealpha',0.1,'edgecolor','none','facecolor','m');
% find the best fit
BL.idx1=inShape(BL.station_bound,stlon1,stlat1);
BL.z1=vertcat(CTD1(BL.idx1).zStar);BL.tr1=vertcat(CTD1(BL.idx1).TR_SF5);
BL.idx2=inShape(BL.station_bound,stlon2,stlat2);
BL.z2=vertcat(CTD2(BL.idx2).zStar);BL.tr2=vertcat(CTD2(BL.idx2).TR_SF5);
BL.cost1=trapz(zHat_grid,([K.c]-nanmean([CTD1(BL.idx1).c_ij],2)).^2);
BL.cost2=trapz(zHat_grid,([K.c]-nanmean([CTD2(BL.idx2).c_ij],2)).^2);
[~,BL.Ridx1]=min(cost1); [~,BL.Ridx2]=min(cost2);

% Define NI shape
NI.lon=[-59.7416-0.2,-58.7766,-60.1020,-60.4160,-60.9392,-60.6601,...
    -59.6253,-59.4858,-59.4044,-59.1719,-60.1485];
NI.lat=[49.6653-0.1,49.5343-0.25,49.2646,49.1029,48.8717,48.7870,...
    48.3401,48.2631,48.5867,49.0412,49.2569];
NI.station_bound=alphaShape(NI.lon',NI.lat',1,'HoleThreshold',15);
[l,ll]=m_ll2xy(NI.lon,NI.lat);
tmpshp=alphaShape(l',ll',0.01);
NI.AS=plot(tmpshp,'facealpha',0.1,'edgecolor','none','facecolor',rgb_x('forest green'));
% find the best fit
NI.idx1=inShape(NI.station_bound,stlon1,stlat1);
NI.z1=vertcat(CTD1(NI.idx1).zStar);NI.tr1=vertcat(CTD1(NI.idx1).TR_SF5);
NI.idx2=inShape(NI.station_bound,stlon2,stlat2);
NI.z2=vertcat(CTD2(NI.idx2).zStar);NI.tr2=vertcat(CTD2(NI.idx2).TR_SF5);
NI.cost1=trapz(zHat_grid,([K.c]-nanmean([CTD1(NI.idx1).c_ij],2)).^2);
NI.cost2=trapz(zHat_grid,([K.c]-nanmean([CTD2(NI.idx2).c_ij],2)).^2);
[~,NI.Ridx1]=min(cost1); [~,NI.Ridx2]=min(cost2);

% Define NA shape
NA.lon=[-59.9509,-60.1369,-60.5904,-60.9392,-61.5205,-62.7181,-63.6250,...
    -63.2297,-61.8461,-61.1019,-60.5671,-60.1834,-60.1020,-60.4160];
NA.lat=[49.5651,49.2800,49.0412,49.1722,49.3648,49.7115,49.9734,50.1044,...
    50.0042,49.9888,50.0505,49.9734,49.2646,49.1029];
NA.station_bound=alphaShape(NA.lon',NA.lat',1,'HoleThreshold',15);
[l,ll]=m_ll2xy(NA.lon,NA.lat);
tmpshp=alphaShape(l',ll',0.01);
NA.AS=plot(tmpshp,'facealpha',0.1,'edgecolor','none','facecolor',rgb_x('turquoise'));
% find the best fit
NA.idx1=inShape(NA.station_bound,stlon1,stlat1);
NA.z1=vertcat(CTD1(NA.idx1).zStar);NA.tr1=vertcat(CTD1(NA.idx1).TR_SF5);
NA.idx2=inShape(NA.station_bound,stlon2,stlat2);
NA.z2=vertcat(CTD2(NA.idx2).zStar);NA.tr2=vertcat(CTD2(NA.idx2).TR_SF5);
NA.cost1=trapz(zHat_grid,([K.c]-nanmean([CTD1(NA.idx1).c_ij],2)).^2);
NA.cost2=trapz(zHat_grid,([K.c]-nanmean([CTD2(NA.idx2).c_ij],2)).^2);
[~,NA.Ridx1]=min(cost1); [~,NA.Ridx2]=min(cost2);

%%%%%%%%%%%%%%%%%%%%% Plot profiles
ax2=axes('color','none','position',[.25 0.725 0.075 0.125 ]);
hold on;
line([0 0.8],[zStar_grid(sigIdx1) zStar_grid(sigIdx1)],...
    'color','k','linestyle','--','linewidth',0.75);
scatter(WL.tr1,WL.z1,5,'filled','markerfacealpha',0.4','markeredgecolor','none',...
    'markerfacecolor','b');
scatter(WL.tr2,WL.z2,5,'filled','markerfacealpha',0.4','markeredgecolor','none',...
    'markerfacecolor','r');

ylim([150 400]);xlim([0 0.8]);xl=xlim;
xticks([0 .25 .5 .75 1]);xtickangle(45);
set(gca,'ydir','reverse');
ylabel('\it{z^*} \rm{(m)}'); xlabel('SF_5CF_3 (fmol kg^{-1})');
grid on; set(gca,'GridColor',[0 0 0],'GridAlpha',0.25);

% Hypsography
idx=inShape(WL.station_bound,bathy.lon,bathy.lat);
ssetZ=bathy.Z(idx); ssetZ(isnan(ssetZ))=[];
ssetA=allA(idx);ssetZ(isnan(ssetA))=[];
hypsZ=0:max(ssetZ); hyps=NaN(length(hypsZ),1);
for i=0:max(hypsZ)
    hyps(i+1,1)=sum(ssetA(ssetZ>i));
end
hyps = (hyps - min(hyps)) / (max(hyps) - min(hyps))* max(xl);
l=area([hyps;0],[hypsZ';max(hypsZ)],'facecolor','none','facealpha',0.2,...
    'edgecolor','k');
uistack(l,'bottom');
r=rectangle('Position', [ax2.XLim(1), ax2.YLim(1), diff(ax2.XLim), diff(ax2.YLim)], 'FaceColor', [1 1 1 0.5], 'EdgeColor', 'none');
uistack(r,'bottom');

%%%%%%%%%%%%%%%%%%%%%
ax3=axes('color','none','position',[.41 0.635 0.075 0.125]); hold on;
line([0 0.8],[zStar_grid(sigIdx1) zStar_grid(sigIdx1)],...
    'color','k','linestyle','--','linewidth',0.75);
scatter(EL.tr1,EL.z1,5,'filled','markerfacealpha',0.4','markeredgecolor','none',...
    'markerfacecolor','b');
scatter(EL.tr2,EL.z2,5,'filled','markerfacealpha',0.4','markeredgecolor','none',...
    'markerfacecolor','r');

ylim([150 400]);xlim([0 0.8]);xl=xlim;
xticks([0 .25 .5 .75 1]);xticklabels([]);yticklabels([]);xtickangle(0);
set(gca,'ydir','reverse');
%ylabel('\it{z^*} \rm{(m)}'); %xlabel('SF5 (fmol kg^{-1})');
grid on; set(gca,'GridColor',[0 0 0],'GridAlpha',0.25);

% Hypsography
idx=inShape(EL.station_bound,bathy.lon,bathy.lat);
ssetZ=bathy.Z(idx); ssetZ(isnan(ssetZ))=[];
ssetA=allA(idx);ssetZ(isnan(ssetA))=[];
hypsZ=0:max(ssetZ); hyps=NaN(length(hypsZ),1);
for i=0:max(hypsZ)
    hyps(i+1,1)=sum(ssetA(ssetZ>i));
end
hyps = (hyps - min(hyps)) / (max(hyps) - min(hyps))* max(xl);
l=area([hyps;0],[hypsZ';max(hypsZ)],'facecolor','none','facealpha',0.2,...
    'edgecolor','k');
uistack(l,'bottom');
r=rectangle('Position', [ax3.XLim(1), ax3.YLim(1), diff(ax3.XLim), diff(ax3.YLim)], 'FaceColor', [1 1 1 0.5], 'EdgeColor', 'none');
uistack(r,'bottom');

%%%%%%%%%%%%%%%%%%%%%
ax4=axes('color','none','position',[.515 0.6 0.075 0.125 ]); hold on;
line([0 0.8],[zStar_grid(sigIdx1) zStar_grid(sigIdx1)],...
    'color','k','linestyle','--','linewidth',0.75);
scatter(OL.tr1,OL.z1,5,'filled','markerfacealpha',0.4','markeredgecolor','none',...
    'markerfacecolor','b');
scatter(OL.tr2,OL.z2,5,'filled','markerfacealpha',0.4','markeredgecolor','none',...
    'markerfacecolor','r');

ylim([150 400]);xlim([0 0.8]);xl=xlim;
xticks([0 .25 .5 .75 1]);xticklabels([]);yticklabels([]);xtickangle(0);
set(gca,'ydir','reverse');
%ylabel('\it{z^*} \rm{(m)}'); %xlabel('SF5 (fmol kg^{-1})');
grid on; set(gca,'GridColor',[0 0 0],'GridAlpha',0.25);

% Hypsography
idx=inShape(OL.station_bound,bathy.lon,bathy.lat);
ssetZ=bathy.Z(idx); ssetZ(isnan(ssetZ))=[];
ssetA=allA(idx);ssetZ(isnan(ssetA))=[];
hypsZ=0:max(ssetZ); hyps=NaN(length(hypsZ),1);
for i=0:max(hypsZ)
    hyps(i+1,1)=sum(ssetA(ssetZ>i));
end
hyps = (hyps - min(hyps)) / (max(hyps) - min(hyps))* max(xl);
l=area([hyps;0],[hypsZ';max(hypsZ)],'facecolor','none','facealpha',0.2,...
    'edgecolor','k');
uistack(l,'bottom');
r=rectangle('Position', [ax4.XLim(1), ax4.YLim(1), diff(ax4.XLim), diff(ax4.YLim)], 'FaceColor', [1 1 1 0.5], 'EdgeColor', 'none');
uistack(r,'bottom');

%%%%%%%%%%%%%%%%%%%%%
ax5=axes('color','none','position',[.75 0.415 0.075 0.125 ]); hold on;
line([0 0.8],[zStar_grid(sigIdx1) zStar_grid(sigIdx1)],...
    'color','k','linestyle','--','linewidth',0.75);
scatter(CS.tr1,CS.z1,5,'filled','markerfacealpha',0.4','markeredgecolor','none',...
    'markerfacecolor','b');
scatter(CS.tr2,CS.z2,5,'filled','markerfacealpha',0.4','markeredgecolor','none',...
    'markerfacecolor','r');

ylim([150 400]);xlim([0 0.8]);xl=xlim;
xticks([0 .25 .5 .75 1]);xticklabels([]);yticklabels([]);xtickangle(0);
set(gca,'ydir','reverse');
%ylabel('\it{z^*} \rm{(m)}'); %xlabel('SF5 (fmol kg^{-1})');
grid on; set(gca,'GridColor',[0 0 0],'GridAlpha',0.25);

% Hypsography
idx=inShape(CS.station_bound,bathy.lon,bathy.lat);
ssetZ=bathy.Z(idx); ssetZ(isnan(ssetZ))=[];
ssetA=allA(idx);ssetZ(isnan(ssetA))=[];
hypsZ=0:max(ssetZ); hyps=NaN(length(hypsZ),1);
for i=0:max(hypsZ)
    hyps(i+1,1)=sum(ssetA(ssetZ>i));
end
hyps = (hyps - min(hyps)) / (max(hyps) - min(hyps))* max(xl);
l=area([hyps;0],[hypsZ';max(hypsZ)],'facecolor','none','facealpha',0.2,...
    'edgecolor','k');
uistack(l,'bottom');
r=rectangle('Position', [ax5.XLim(1), ax5.YLim(1), diff(ax5.XLim), diff(ax5.YLim)], 'FaceColor', [1 1 1 0.5], 'EdgeColor', 'none');
uistack(r,'bottom');

%%%%%%%%%%%%%%%%%%%%%
ax7=axes('color','none','position',[.75 0.775 0.075 0.125 ]);hold on;
line([0 0.8],[zStar_grid(sigIdx1) zStar_grid(sigIdx1)],...
    'color','k','linestyle','--','linewidth',0.75);
scatter(BL.tr1,BL.z1,5,'filled','markerfacealpha',0.4','markeredgecolor','none',...
    'markerfacecolor','b');
scatter(BL.tr2,BL.z2,5,'filled','markerfacealpha',0.4','markeredgecolor','none',...
    'markerfacecolor','r');

ylim([150 400]);xlim([0 0.8]);xl=xlim;
xticks([0 .25 .5 .75 1]);xticklabels([]);yticklabels([]);xtickangle(0);
set(gca,'ydir','reverse');
%ylabel('\it{z^*} \rm{(m)}'); %xlabel('SF5 (fmol kg^{-1})');
grid on; set(gca,'GridColor',[0 0 0],'GridAlpha',0.25);

% Hypsography
idx=inShape(BL.station_bound,bathy.lon,bathy.lat);
ssetZ=bathy.Z(idx); ssetZ(isnan(ssetZ))=[];
ssetA=allA(idx);ssetZ(isnan(ssetA))=[];
hypsZ=0:max(ssetZ); hyps=NaN(length(hypsZ),1);
for i=0:max(hypsZ)
    hyps(i+1,1)=sum(ssetA(ssetZ>i));
end
hyps = (hyps - min(hyps)) / (max(hyps) - min(hyps))* max(xl);
l=area([hyps;0],[hypsZ';max(hypsZ)],'facecolor','none','facealpha',0.2,...
    'edgecolor','k');
uistack(l,'bottom');
r=rectangle('Position', [ax7.XLim(1), ax7.YLim(1), diff(ax7.XLim), diff(ax7.YLim)], 'FaceColor', [1 1 1 0.5], 'EdgeColor', 'none');
uistack(r,'bottom');

%%%%%%%%%%%%%%%%%%%%%
ax8=axes('color','none','position',[.675 0.65 0.065 0.115 ]); hold on;
line([0 0.8],[zStar_grid(sigIdx1) zStar_grid(sigIdx1)],...
    'color','k','linestyle','--','linewidth',0.75);
scatter(NI.tr1,NI.z1,5,'filled','markerfacealpha',0.4','markeredgecolor','none',...
    'markerfacecolor','b');
scatter(NI.tr2,NI.z2,5,'filled','markerfacealpha',0.4','markeredgecolor','none',...
    'markerfacecolor','r');

ylim([150 400]);xlim([0 0.8]);xl=xlim;
xticks([0 .25 .5 .75 1]);xticklabels([]);yticklabels([]);xtickangle(0);
set(gca,'ydir','reverse');
%ylabel('\it{z^*} \rm{(m)}'); %xlabel('SF5 (fmol kg^{-1})');
grid on; set(gca,'GridColor',[0 0 0],'GridAlpha',0.25);

% Hypsography
idx=inShape(NI.station_bound,bathy.lon,bathy.lat);
ssetZ=bathy.Z(idx); ssetZ(isnan(ssetZ))=[];
ssetA=allA(idx);ssetZ(isnan(ssetA))=[];
hypsZ=0:max(ssetZ); hyps=NaN(length(hypsZ),1);
for i=0:max(hypsZ)
    hyps(i+1,1)=sum(ssetA(ssetZ>i));
end
hyps = (hyps - min(hyps)) / (max(hyps) - min(hyps))* max(xl);
l=area([hyps;0],[hypsZ';max(hypsZ)],'facecolor','none','facealpha',0.2,...
    'edgecolor','k');
uistack(l,'bottom');
r=rectangle('Position', [ax8.XLim(1), ax8.YLim(1), diff(ax8.XLim), diff(ax8.YLim)], 'FaceColor', [1 1 1 0.5], 'EdgeColor', 'none');
uistack(r,'bottom');

%%%%%%%%%%%%%%%%%%%%%
ax9=axes('color','none','position',[.55 0.75 0.075 0.125 ]); hold on;
line([0 0.8],[zStar_grid(sigIdx1) zStar_grid(sigIdx1)],...
    'color','k','linestyle','--','linewidth',0.75);
scatter(NA.tr1,NA.z1,5,'filled','markerfacealpha',0.4','markeredgecolor','none',...
    'markerfacecolor','b');
scatter(NA.tr2,NA.z2,5,'filled','markerfacealpha',0.4','markeredgecolor','none',...
    'markerfacecolor','r');

ylim([150 400]);xlim([0 0.8]);xl=xlim;
xticks([0 .25 .5 .75 1]);xticklabels([]);yticklabels([]);xtickangle(0);
set(gca,'ydir','reverse');
%ylabel('\it{z^*} \rm{(m)}'); %xlabel('SF5 (fmol kg^{-1})');
grid on; set(gca,'GridColor',[0 0 0],'GridAlpha',0.25);

% Hypsography
idx=inShape(NA.station_bound,bathy.lon,bathy.lat);
ssetZ=bathy.Z(idx); ssetZ(isnan(ssetZ))=[];
ssetA=allA(idx);ssetZ(isnan(ssetA))=[];
hypsZ=0:max(ssetZ); hyps=NaN(length(hypsZ),1);
for i=0:max(hypsZ)
    hyps(i+1,1)=sum(ssetA(ssetZ>i));
end
hyps = (hyps - min(hyps)) / (max(hyps) - min(hyps))* max(xl);
l=area([hyps;0],[hypsZ';max(hypsZ)],'facecolor','none','facealpha',0.2,...
    'edgecolor','k');
uistack(l,'bottom');
r=rectangle('Position', [ax9.XLim(1), ax9.YLim(1), diff(ax9.XLim), diff(ax9.YLim)], 'FaceColor', [1 1 1 0.5], 'EdgeColor', 'none');
uistack(r,'bottom');


%%%%%%%%%%%%%%%
iax=axes('color','w','position',[0.165 0.415 0.15+0.165 0.2]); hold on;
grid on; axis ij;ylim([-125 125]); box on; xticks([]);

axV1=axes('color','w','position',[0.165 0.415 0.15 0.2 ]); hold on;
%%%% CIOPS-E
plot(movmean(model10.Cp/sum(model10.Cp),10),model10.zco,'color',...
    rgb_x('light blue'),'linestyle','-','linewidth',1);
ll=plot(K(CidxM).rM1,zHat_grid,'-.','color',rgb_x('light blue'),'linewidth',1);
Color=[ll.Color 0.1];

%%%% Leg 1
l(1)=plot(K(Cidx1).rM1,zHat_grid,'color',rgb_x('light blue'),'linewidth',2);
l(1).Color=[l(1).Color 0.7];
plot(ri1,zHat_grid,'b--','linewidth',1);
errorbar(ri1(1:25:length(zHat_grid)),zHat_grid(1:25:length(zHat_grid)),...
    errLim1(1:25:length(zHat_grid)),...
    '.','horizontal','markersize',5,'marker','none',...
    'markerfacecolor','b','markeredgecolor','none',...
    'color','b','CapSize',2,'linewidth',0.5);
scatter(ri1(1:25:length(zHat_grid)),zHat_grid(1:25:length(zHat_grid)),8,'s',...
    'markerfacecolor','b','markeredgecolor','b','markerfacealpha',0.7);

%%%% Leg 2
l(2)=plot(K(Cidx2).rM2,zHat_grid,'-','color',rgb_x('rose'),'linewidth',2);
l(2).Color=[l(2).Color 0.7];
plot(ri2,zHat_grid,'r--','linewidth',1);
errorbar(ri2(10:25:length(zHat_grid)),zHat_grid(10:25:length(zHat_grid)),...
    errLim2(10:25:length(zHat_grid))*2,...
    '.','horizontal','markersize',5,'marker','none',...
    'markerfacecolor','r','markeredgecolor','none',...
    'color','r','CapSize',2,'linewidth',0.5);
scatter(ri2(10:25:length(zHat_grid)),zHat_grid(10:25:length(zHat_grid)),8,'s',...
    'markerfacecolor','r','markeredgecolor','r','markerfacealpha',0.7);
xlabel('\it{r}');
ylabel('$$\hat{z}\, (m)$$', 'Interpreter', 'latex');
set(gcf,'color','w')
grid on; axis tight; axis ij;ylim([-125 125]);
text(0.021,-110,'b)','fontsize',8);

%%%% Depth segregation %%%%%%
axV2=axes('color','w','position',[0.165+0.165 0.415 0.15 0.2]); hold on;
l=plot(K(CidxDeep1).rM1,zHat_grid,'linewidth',2);
l.Color=[rgb_x('cerulean') 0.5];
l=plot(K(CidxShallow1).rM1,zHat_grid,'linewidth',2);
l.Color=[rgb_x('algae green') 0.5];

l=plot(K(CidxDeep2).rM2,zHat_grid,'linewidth',2);
l.Color=[rgb_x('cerulean').*0.8 0.7];
l=plot(K(CidxShallow2).rM2,zHat_grid,'linewidth',2);
l.Color=[rgb_x('algae green').*0.8 0.7];

plot(smooth(deep_ri1,10),zHat_grid,'--','color',rgb_x('cerulean'));
errorbar(deep_ri1(10:25:length(zHat_grid)),zHat_grid(10:25:length(zHat_grid)),...
    errDeep1(10:25:length(zHat_grid)),...
    '.','horizontal','markersize',5,'marker','none',...
    'markerfacecolor',rgb_x('cerulean'),'markeredgecolor','none',...
    'color',rgb_x('cerulean'),'CapSize',2,'linewidth',0.5);
scatter(deep_ri1(10:25:length(zHat_grid)),zHat_grid(10:25:length(zHat_grid)),8,'s',...
    'markerfacecolor',rgb_x('cerulean'),'markeredgecolor',rgb_x('cerulean'),'markerfacealpha',0.5);

plot(smooth(shallow_ri1,10),zHat_grid,'--','color',rgb_x('algae green'));
errorbar(shallow_ri1(10:25:length(zHat_grid)),zHat_grid(10:25:length(zHat_grid)),...
    errShallow1(10:25:length(zHat_grid)),...
    '.','horizontal','markersize',5,'marker','none',...
    'markerfacecolor',rgb_x('algae green'),'markeredgecolor','none',...
    'color',rgb_x('algae green'),'CapSize',2,'linewidth',0.5);
scatter(shallow_ri1(10:25:length(zHat_grid)),zHat_grid(10:25:length(zHat_grid)),8,'s',...
    'markerfacecolor',rgb_x('algae green'),'markeredgecolor',rgb_x('algae green'),'markerfacealpha',0.5);

plot(smooth(deep_ri2,10),zHat_grid,'--','color',rgb_x('cerulean'));
errorbar(deep_ri2(10:25:length(zHat_grid)),zHat_grid(10:25:length(zHat_grid)),...
    errDeep2(10:25:length(zHat_grid)),...
    '.','horizontal','markersize',5,'marker','none',...
    'markerfacecolor',rgb_x('cerulean'),'markeredgecolor','none',...
    'color',rgb_x('cerulean'),'CapSize',2,'linewidth',0.5);
scatter(deep_ri2(10:25:length(zHat_grid)),zHat_grid(10:25:length(zHat_grid)),8,'s',...
    'markerfacecolor',rgb_x('cerulean').*0.8,'markeredgecolor',rgb_x('cerulean'),'markerfacealpha',0.5);

plot(smooth(shallow_ri2,10),zHat_grid,'--','color',rgb_x('algae green'));
errorbar(shallow_ri2(10:25:length(zHat_grid)),zHat_grid(10:25:length(zHat_grid)),...
    errShallow2(10:25:length(zHat_grid)),...
    '.','horizontal','markersize',5,'marker','none',...
    'markerfacecolor',rgb_x('algae green'),'markeredgecolor','none',...
    'color',rgb_x('algae green'),'CapSize',2,'linewidth',0.5);
scatter(shallow_ri2(10:25:length(zHat_grid)),zHat_grid(10:25:length(zHat_grid)),8,'s',...
    'markerfacecolor',rgb_x('algae green').*0.8,'markeredgecolor',rgb_x('algae green'),'markerfacealpha',0.5);

axis ij
set(gcf,'color','w');
grid on; axis tight; axis ij;ylim([-125 125]);yticklabels([]);
axV2.GridColor=[0.15 0.15 0.15];
axV2.YColor=[0.7 0.7 0.7];
text(0.021,-110,'c)','fontsize',8);

set(findall(gcf,'-property','fontname'),'fontname','Latin Modern Roman');
set(findall(gcf,'-property','fontsize'),'fontsize',8);


%%
export_fig /Users/samst/Dropbox/UBC/TReX_paper/figures/letter/verticalV2.pdf -dpdf -nofontswap

%% Scaling analysis (Cyr equation 5)
% Not CYR e5
shallowpp=sum(bathy.Z>150 & bathy.Z<300,1:2)
deeppp=sum(bathy.Z>300,1:2)


