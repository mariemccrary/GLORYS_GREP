% Marie McCrary April 20, 2020

% This file loads temperature data from monthly GLORYS GREP v2, subsets it
% by polygon mask, and calculates heat content within that mask.

% The first data file contans the data for each of the reanalyses. The
% second data file contains the mean and standard deviation of variables.

% Several calculations include functions from the Climate Data Toolbox:
% Chad A. Greene, Kaustubh Thirumalai, Kelly A. Kearney, José Miguel
% Delgado, Wolfgang Schwanghart, Natalie S. Wolfenbarger, Kristen M. Thyng,
% David E. Gwyther, Alex S. Gardner, and Donald D. Blankenship. The Climate
% Data Toolbox for MATLAB. Geochemistry, Geophysics, Geosystems 2019.
% doi:10.1029/2019GC008392 https://doi.org/10.1029/2019GC008392

data = ('global-reanalysis-phy-001-031-grepv2-monthly_temp.nc');
data_mnstd = ('global-reanalysis-phy-001-031-grepv2-mnstd-monthly.nc');

c = 4186; % J/kgC
rho = 1027; % kg/m3

% Full grid space
lats = double(ncread(data, 'latitude'));
lons = double(ncread(data, 'longitude'));
depth = double(ncread(data, 'depth'));

[lat, lon] = meshgrid(lats, lons);

% Temperature data (longitude,latitude,depth,time)

t_oras = ncread(data, 'thetao_oras');
t_cglo = ncread(data, 'thetao_cglo');
t_glor = ncread(data, 'thetao_glor');
t_foam = ncread(data, 'thetao_foam');

t_mean = ncread(data_mnstd, 'thetao_mean'); 

t_std = ncread(data_mnstd, 'thetao_std');

% Time variable
time = ncread(data, 'time'); % days since 1950-01-01
timespan = length(time);
base = datenum(1950,01,01);
time = (time + base);
time = datenum(time);

date = datevec(time);

%% Grid cell volume

A = cdtarea(lat, lon);

Vol = zeros(length(lons), length(lats), length(depth));
Vol_temp = zeros(length(lons), length(lats), length(depth));

for i = 1:length(depth)
    Vol(:,:,i) = A .* depth(i);
end

for i = 1:8
    Vol_temp(:,:,i) = Vol(:,:,i);
end
Vol10 = nansum(Vol_temp,3);
clear Vol_temp

%% Depth Avg Temperature

d = zeros(length(depth),1);

d(1) = 2*depth(1);
dtot = sum(d);

for i=2:25
     d(i) = 2*(depth(i)-dtot);
     dtot = sum(d);
end

Tsurf(:,:,1,:) = t_mean(:,:,1,:);
Tsurf = squeeze(Tsurf);

% Depth Avg Temperature to 9.8 m, z=8

T10_mean = zeros(length(lons), length(lats), length(depth), timespan);
for i=1:8
    T10_mean(:,:,i,:) = (t_mean(:,:,i,:).*d(i));
end
 
T10d = nansum(T10_mean,3); % Sums T along 3rd dim (depth)
T10d = squeeze(T10d);
Tbar10_mean = T10d./depth(8);

%% Heat content in each layer

% Heat is defined as Q = m*c*deltaT = rho*Vol*c*deltaT

Q_surf = zeros(length(lons), length(lats), timespan);
for i = 1:timespan
    Q_surf(:,:,i) = rho.*c.* Vol(:,:,1) .*Tsurf(:,:,i);
end

Q_10 = zeros(length(lons), length(lats), timespan);

for i = 1:timespan
    Q_10(:,:,i) = rho.*c.* Vol10(:,:) .*Tbar10_mean(:,:,i);
end


%% Subset data by location

% Mackenzie
val = [-136.3858148812 69.0509501902;  -137.9436868227 69.3128192368 ;  -139.3367765793  69.6737303448 ;  -139.4190164954  70.4112187033 ;  -138.9257281747  71.0003029587 ;  -137.0775977055 71.4307629332 ;  -135.7256279071  71.3044102369 ;  -132.8894090333  71.0957989552 ;  -132.3757119098 70.4830637218 ;  -131.9398705894 69.8884664616 ;  -132.9863130522 69.6935711514 ;  -133.8059910402 69.5899640842 ;  -134.4672362486 69.7206254855 ;  -134.915474027 69.5695514246 ;  -135.304753042 69.4292489187 ;  -135.8775601052 69.2859396466 ;-136.3858148812 69.0509501902];
X = val(:, 2);
Y = val(:, 1);

bound = geoshape(X,Y);
mask = geomask(lat, lon, bound.Latitude, bound.Longitude);

Tbar10_M = zeros(length(lons), length(lats), timespan);

for i=1:timespan
    Tbar10_M(:,:,i) = Tbar10_mean(:,:,i).*mask(:,:);
end

Q_10_M = zeros(length(lons), length(lats), timespan);

for i = 1:timespan
    Q_10_M(:,:,i) = rho.*c.* Vol10(:,:) .*Tbar10_M(:,:,i);
end

Q_10_M_sum = sum(Q_10_M, [1 2]);
Q_10_M_sum = squeeze(Q_10_M_sum);

figure
hold on
plot(time, Q_10_M_sum)
dateFormat = 11;
datetick('x','mm/yy')

%% Plot full data and subsetted data to compare

figure
hold on
worldmap([65 80], [-180 -130])
pcolorm(lat, lon, (Q_10_M(:,:,57))./(10^17));
cmocean('thermal')
caxis([0 6])
h=colorbar;
h.Location = 'southoutside';
geoshow('landareas.shp', 'FaceColor', [0.85 0.85 0.85])
geoshow('worldrivers.shp','Color', 'blue')
ylabel(h,'J (x 10^{17})');
mlabel('off');
plabel('off');
title('10 m','FontWeight','normal')
%% NOTES %%
%Tbar, Htotal and 
%estimate, FWtotal, Sbar, Rhobar
% h=colorbar('SouthOutside');
% set(h, 'Position', [.1 .05 .8150 .05]);
% h.Label.String= 'Sea Ice Concentration (%)';
% %% Ice Data
% [ci,lat,lon] = arcticseaice('September 15, 2007', 'km', 'noplot');
% 
% figure
% hold on
% worldmap([65 90], [-180 -130])
% pcolorm(lat, lon, ci)
% cmocean('ice')
% h=colorbar;
% caxis([0 100])
% h.Location = 'southoutside';
% geoshow('landareas.shp', 'FaceColor', [0.85 0.85 0.85])
% geoshow('worldrivers.shp','Color', 'blue')
% hold on
% ylabel(h,'Concentration (%)')
% mlabel('off');
% plabel('off');
% title('Ice','FontWeight','normal')
% 
% % savefig('')

