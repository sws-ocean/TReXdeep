%% Extract east coast bathymetry
clear
addpath(genpath('/Users/samst/Dropbox/UBC/TReX_paper/'));

lat_lim=[46.75 51.5];
lon_lim=[-69.5 -57.5];

allLon=ncread('eastern_canada_3_isl_2015.nc','x');
allLat=ncread('eastern_canada_3_isl_2015.nc','y');

lonIdx=find(allLon>lon_lim(1) & allLon<lon_lim(2));
latIdx=find(allLat>lat_lim(1) & allLat<lat_lim(2));

% Load in decimated data
Z=ncread('eastern_canada_3_isl_2015.nc','z',...
    [min(lonIdx) min(latIdx)],[length(lonIdx) length(latIdx)]);

% Decimate
dec=2;
tmplon=allLon(lonIdx((1:dec:end)));tmplat=allLat(latIdx((1:dec:end)));
[lat,lon]=meshgrid(tmplat,tmplon);
Z=Z(1:dec:end,1:dec:end);
Z(isnan(Z))=NaN;

figure
contourf(lon,lat,Z,0:-100:-400,'linestyle','none');
 
save('/Users/samst/Dropbox/UBC/TReX_paper/data/GoSL_Zhires.mat',...
    'lon','lat','Z');
